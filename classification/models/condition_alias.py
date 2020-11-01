import operator
from functools import reduce
from typing import List, Set, Optional

from django.db import models, transaction
from django.db.models import Count, Q
from django.db.models.signals import post_save
from django.dispatch import receiver
from django_extensions.db.models import TimeStampedModel
from django.contrib.auth.models import User
from guardian.shortcuts import assign_perm
from classification.regexes import db_ref_regexes
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from snpdb.models import Lab


class ConditionAliasStatus(models.TextChoices):
    PENDING = 'P', 'Pending'
    RESOLVED_AUTO = 'A', 'Resolved (Auto)'
    RESOLVED = 'R', 'Resolved'
    UNMATCHABLE = 'U', 'Unmatchable'


class ConditionAliasJoin(models.TextChoices):
    OR = 'O', 'Or'  # aka uncertain
    AND = 'A', 'And'  # aka combined


class ConditionAlias(TimeStampedModel, GuardianPermissionsMixin):
    lab = models.ForeignKey(to=Lab, on_delete=models.CASCADE)
    source_text = models.TextField()
    source_gene_symbol = models.TextField()
    status = models.CharField(max_length=1, choices=ConditionAliasStatus.choices, default=ConditionAliasStatus.PENDING)
    join_mode = models.CharField(max_length=1, choices=ConditionAliasJoin.choices, default=ConditionAliasJoin.OR)
    records_affected = models.IntegerField()
    aliases = models.JSONField(null=True, blank=True)
    updated_by = models.ForeignKey(to=User, on_delete=models.SET_NULL, null=True, blank=True)

    class Meta:
        unique_together = ("lab", "source_text", "source_gene_symbol")

    def __str__(self):
        alias_text = "?"
        if self.aliases:
            alias_text = f" {self.get_join_mode_display().lower()} ".join(self.aliases)
        return f"{self.lab.name} | {self.source_gene_symbol} | {self.source_text} -> {alias_text}"

    @staticmethod
    def sync_aliases():
        from classification.models import ClassificationModification
        qs = ClassificationModification.objects.filter(is_last_published=True)\
                .filter(published_evidence__condition__db_refs__isnull=True)\
                .values("classification__lab", "published_evidence__condition__value", "published_evidence__gene_symbol__value")\
                .annotate(total=Count('published_evidence__condition__value'))\
                .order_by("classification__lab", "published_evidence__condition__value", "published_evidence__gene_symbol__value")

        ConditionAlias.objects.update(records_affected=0)
        for result in qs:
            if condition := result.get('published_evidence__condition__value'):
                if gene_symbol := result.get('published_evidence__gene_symbol__value'):
                    lab_id = result.get('classification__lab')
                    # double check that we can't find some terms automatically
                    if not db_ref_regexes.search(text=condition):
                        total = result.get('total')

                        record, _ = ConditionAlias.objects.update_or_create(
                            lab_id = lab_id,
                            source_text = condition,
                            source_gene_symbol = gene_symbol,
                            defaults={
                                "records_affected": total
                            }
                        )
                        record.attempt_auto_match()

    @staticmethod
    def _name_to_words(name: str) -> Set[str]:
        cleaned = name.replace(",", " ").lower()
        words = cleaned.split(" ")
        words = [word.strip() for word in words]
        return set([word for word in words if word])

    def attempt_auto_match(self):
        from annotation.models import MonarchDiseaseOntology
        word_set = ConditionAlias._name_to_words(self.source_text)
        q_list: List[Q] = [Q(name__icontains=word) for word in word_set]
        mondos = MonarchDiseaseOntology.objects.filter(reduce(operator.and_, q_list))

        best_match: Optional[MonarchDiseaseOntology] = None

        mondo: MonarchDiseaseOntology
        for mondo in mondos:
            defn = mondo.definition
            if defn:
                defn = defn.lower()
            leftover_words = ConditionAlias._name_to_words(mondo.name)
            for term in word_set:
                if term in leftover_words:  # it should always be a term due to the query
                    leftover_words.remove(term)
            if not leftover_words:
                best_match = term
                break
            if len(leftover_words) == 1 and self.source_gene_symbol in defn:
                leftover_word = list(leftover_words)[0]
                # want to see if single leftover word is number?
                print(f"Leftover word = {leftover_word}")
                best_match = term
                break

        if best_match:
            self.aliases = [best_match.id_str]
            self.status = ConditionAliasStatus.RESOLVED_AUTO
            self.save()
            return True
        else:
            return False


@receiver(post_save, sender=ConditionAlias)
def set_condition_alias_permissions(sender, created: bool, instance: ConditionAlias, **kwargs):
    if created:
        group = instance.lab.group
        assign_perm(ConditionAlias.get_read_perm(), group, instance)
        assign_perm(ConditionAlias.get_write_perm(), group, instance)


class ConditionAliasSearchCache(TimeStampedModel):
    search_text = models.TextField(primary_key=True)
    cached_results = models.JSONField()
