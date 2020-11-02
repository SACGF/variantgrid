import operator
from dataclasses import dataclass
from functools import reduce
from typing import List, Set, Optional

from django.db import models
from django.db.models import Count, Q
from django.db.models.signals import post_save
from django.dispatch import receiver
from django_extensions.db.models import TimeStampedModel
from django.contrib.auth.models import User
from guardian.shortcuts import assign_perm

from annotation.models import MonarchDiseaseOntology
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


@dataclass(frozen=True)
class ConditionAliasAutoMatch:
    contains_all_count: int
    best_match: Optional[MonarchDiseaseOntology]
    secondary_matches: List[MonarchDiseaseOntology]

    def __str__(self):
        secondary_str = self.secondary_matches[0].name if len(self.secondary_matches) == 1 else str(len(self.secondary_matches))
        did_match = self.best_match is not None or len(self.secondary_matches) == 1
        return f"Matched = {did_match} - {self.contains_all_count} candidate, - perfect match = {self.best_match.name if self.best_match else 'None'} secondary = {secondary_str}"


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
                            lab_id=lab_id,
                            source_text=condition,
                            source_gene_symbol=gene_symbol,
                            defaults={
                                "records_affected": total
                            }
                        )
                        match_summary = record.attempt_auto_match()
                        print(f"{record.source_text} -> {match_summary}")

    @staticmethod
    def _name_to_words(name: str) -> Set[str]:
        cleaned = name.replace(",", " ").lower()
        words = cleaned.split(" ")
        words = [word.strip() for word in words]
        return set([word for word in words if word])

    def attempt_auto_match(self) -> ConditionAliasAutoMatch:
        word_set = ConditionAlias._name_to_words(self.source_text)
        q_list: List[Q] = [Q(name__icontains=word) for word in word_set]
        mondos = list(MonarchDiseaseOntology.objects.filter(reduce(operator.and_, q_list)).order_by('name'))

        possible_matches_count = len(mondos)

        best_match: Optional[MonarchDiseaseOntology] = None
        secondary_matches: List[MonarchDiseaseOntology] = list()
        gene_symbol = self.source_gene_symbol.lower() if self.source_gene_symbol else None

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
                best_match = mondo
            elif len(leftover_words) <= 2 and gene_symbol and defn and gene_symbol in defn:
                secondary_matches.append(mondo)

        use_match: Optional[MonarchDiseaseOntology]
        if best_match:
            use_match = best_match
        elif len(secondary_matches) == 1:
            use_match = secondary_matches[0]

        if best_match:
            self.aliases = [best_match.id_str]
            self.status = ConditionAliasStatus.RESOLVED_AUTO
            self.save()

        return ConditionAliasAutoMatch(
            contains_all_count=possible_matches_count,
            best_match=best_match,
            secondary_matches=secondary_matches
        )


@receiver(post_save, sender=ConditionAlias)
def set_condition_alias_permissions(sender, created: bool, instance: ConditionAlias, **kwargs):
    if created:
        group = instance.lab.group
        assign_perm(ConditionAlias.get_read_perm(), group, instance)
        assign_perm(ConditionAlias.get_write_perm(), group, instance)


class ConditionAliasSearchCache(TimeStampedModel):
    search_text = models.TextField(primary_key=True)
    cached_results = models.JSONField()
