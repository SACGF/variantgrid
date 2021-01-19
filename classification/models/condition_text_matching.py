from typing import List, Optional, Dict, Any, Tuple, Set
import re

from django.contrib.auth.models import User
from django.contrib.postgres.fields import ArrayField
from django.db.models import CASCADE, QuerySet, SET_NULL
from django.dispatch import receiver
from guardian.shortcuts import assign_perm
from lazy import lazy
from model_utils.models import TimeStampedModel
from django.db import models
from classification.enums import SpecialEKeys, ShareLevel
from classification.models import Classification, ClassificationModification, classification_post_publish_signal, \
    flag_types
from classification.regexes import db_ref_regexes, DbRegexes
from flags.models import flag_comment_action, Flag, FlagComment, FlagResolution
from genes.models import GeneSymbol
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from ontology.ontology_matching import OntologyMatching, normalize_condition_text
from snpdb.models import Lab


class ConditionText(TimeStampedModel, GuardianPermissionsMixin):
    normalized_text = models.TextField()
    lab = models.ForeignKey(Lab, on_delete=CASCADE, null=True, blank=True)

    # set to true if root condition text was set by auto-matching
    requires_approval = models.BooleanField(blank=True, default=False)
    last_edited_by = models.ForeignKey(User, on_delete=SET_NULL, null=True, blank=True)

    classifications_count = models.IntegerField(default=0)
    classifications_count_outstanding = models.IntegerField(default=0)

    class Meta:
        unique_together = ("normalized_text", "lab")

    @staticmethod
    def normalize(text: str):
        return normalize_condition_text(text)

    def __str__(self):
        return self.normalized_text

    # TODO is this better done as an before save hook?
    def save(self, **kwargs):

        super().save(**kwargs)
        assign_perm(self.get_read_perm(), self.lab.group, self)
        assign_perm(self.get_write_perm(), self.lab.group, self)


class MultiCondition(models.TextChoices):
    NOT_DECIDED = 'X', 'Not decided'
    UNCERTAIN = 'U', 'Uncertain'  # aka uncertain
    CO_OCCURRING = 'C', 'Co-occurring'  # aka combined


class ResolvedCondition:

    def __init__(self):
        self.condition_multi_operation = MultiCondition.NOT_DECIDED
        self.condition_xrefs = OntologyMatching()


class ConditionTextMatch(TimeStampedModel, GuardianPermissionsMixin):
    condition_text = models.ForeignKey(ConditionText, on_delete=CASCADE)

    parent = models.ForeignKey('ConditionTextMatch', on_delete=CASCADE, null=True, blank=True)
    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE, null=True, blank=True)
    mode_of_inheritance = ArrayField(models.TextField(blank=False), default=None, null=True, blank=True)
    classification = models.OneToOneField(Classification, on_delete=CASCADE, null=True, blank=True)

    @classmethod
    def get_permission_object(self):
        return self.condition_text

    @classmethod
    def get_permission_class(cls):
        return ConditionText

    @property
    def hierarchy(self) -> List[Any]:
        h_list = list()
        if self.gene_symbol:
            h_list.append(self.gene_symbol)
        if self.mode_of_inheritance is not None:
            h_list.append(self.mode_of_inheritance)
        if self.classification:
            h_list.append(self.classification)
        return h_list

    @property
    def leaf(self) -> Optional[Any]:
        if self.classification:
            return self.classification
        elif self.gene_symbol:
            return self.gene_symbol
        elif self.mode_of_inheritance is not None:
            return self.mode_of_inheritance
        else:
            return None

    # resolved to this,
    condition_xrefs = ArrayField(models.TextField(blank=False), default=list)
    condition_multi_operation = models.CharField(max_length=1, choices=MultiCondition.choices, default=MultiCondition.NOT_DECIDED)

    @lazy
    def resolve_condition_xrefs(self) -> Optional[ResolvedCondition]:
        if not self.condition_xrefs:
            if not self.parent:
                return None
            else:
                return self.parent.resolve_condition_xrefs

        rc = ResolvedCondition()
        rc.condition_multi_operation = self.condition_multi_operation
        for term in self.condition_xrefs:
            rc.condition_xrefs.select_term(term)
        return rc

    @property
    def condition_matching_str(self) -> str:
        result = ''
        if self.condition_xrefs:
            result = ", ".join(self.condition_xrefs)
        if len(self.condition_xrefs) >= 2:
            if self.condition_multi_operation == MultiCondition.NOT_DECIDED:
                result += "; uncertain/co-occurring"
            elif self.condition_multi_operation == MultiCondition.CO_OCCURRING:
                result += "; co-occurring"
            elif self.condition_multi_operation == MultiCondition.UNCERTAIN:
                result += "; uncertain"
        return result

    def update_with_condition_matching_str(self, text: str):
        terms = list()
        condition_multi_operation = MultiCondition.NOT_DECIDED

        parts = text.split(";")
        if len(parts) >= 1:
            terms_part = parts[0]
            # this way we only recognise ids we recognise
            # but do we want that?
            # also allows us to fix the length of ids

            # also assume any dangling number is a mondo term?
            db_refs = db_ref_regexes.search(terms_part, default_regex=DbRegexes.MONDO)
            terms = [db_ref.id_fixed for db_ref in db_refs]
        if len(parts) >= 2:
            operation_part = parts[1].lower()
            if '/' in operation_part:
                condition_multi_operation = MultiCondition.NOT_DECIDED
            elif 'co' in operation_part:
                condition_multi_operation = MultiCondition.CO_OCCURRING
            elif 'un' in operation_part:
                condition_multi_operation = MultiCondition.UNCERTAIN

        self.condition_xrefs = terms
        self.condition_multi_operation = condition_multi_operation

    class Meta:
        unique_together = ("condition_text", "gene_symbol", "mode_of_inheritance", "classification")

    @property
    def is_valid(self):
        return len(self.condition_xrefs) == 1 or \
               (len(self.condition_xrefs) > 1 and self.condition_multi_operation != MultiCondition.NOT_DECIDED)

    @property
    def is_blank(self):
        return len(self.condition_xrefs) == 0

    @lazy
    def children(self) -> QuerySet:
        return self.conditiontextmatch_set.all()

    @staticmethod
    def sync_all():
        cms = ClassificationModification.objects.filter(is_last_published=True, classification__withdrawn=False).select_related("classification", "classification__lab")
        cm: ClassificationModification
        for cm in cms:
            ConditionTextMatch.sync_condition_text_classification(cm=cm, update_counts=False)

        for ct in ConditionText.objects.all():
            update_condition_text_match_counts(ct)
            ct.save()

    @staticmethod
    def sync_condition_text_classification(cm: ClassificationModification, update_counts=True):
        classification = cm.classification

        if classification.withdrawn:
            ConditionTextMatch.objects.filter(classification=classification).delete()
            return

        lab = classification.lab
        gene_str = cm.get(SpecialEKeys.GENE_SYMBOL)
        gene_symbol = GeneSymbol.objects.filter(symbol=gene_str).first()

        existing: ConditionTextMatch = ConditionTextMatch.objects.filter(classification=classification).first()

        if not gene_symbol or classification.withdrawn:
            if existing:
                ct = existing.condition_text
                existing.delete()
                if update_counts:
                    update_condition_text_match_counts(ct)
            return
        else:
            raw_condition_text = cm.get(SpecialEKeys.CONDITION) or ""
            normalized = ConditionText.normalize(raw_condition_text)

            # need mode_of_inheritance to be not null
            mode_of_inheritance = cm.get(SpecialEKeys.MODE_OF_INHERITANCE) or []

            ct: ConditionText
            ct, ct_is_new = ConditionText.objects.get_or_create(normalized_text=normalized, lab=lab)

            # if condition text has changed, remove the old entries
            ConditionTextMatch.objects.filter(classification=classification).exclude(condition_text=ct).delete()

            matches = db_ref_regexes.search(text=normalized)
            root, _ = ConditionTextMatch.objects.get_or_create(
                condition_text=ct,
                gene_symbol=None,
                mode_of_inheritance=None,
                classification=None,
                defaults={"condition_xrefs": [match.id_fixed for match in matches]}
            )

            gene_level, _ = ConditionTextMatch.objects.get_or_create(
                parent=root,
                condition_text=ct,
                gene_symbol=gene_symbol,
                mode_of_inheritance=None,
                classification=None
            )

            mode_of_inheritance_level, _ = ConditionTextMatch.objects.get_or_create(
                parent=gene_level,
                condition_text=ct,
                gene_symbol=gene_symbol,
                mode_of_inheritance=mode_of_inheritance,
                classification=None
            )

            if existing:
                if existing.parent != mode_of_inheritance_level or \
                        existing.condition_text != ct or \
                        existing.gene_symbol != gene_symbol or \
                        existing.mode_of_inheritance != mode_of_inheritance:

                    # update existing to new hierarchy
                    # assume if a condition has been set for this classification specifically that it's
                    # still valid
                    old_root = existing.condition_text
                    existing.parent = mode_of_inheritance_level
                    existing.condition_text = ct
                    existing.mode_of_inheritance = mode_of_inheritance
                    existing.save()

                    if update_counts:
                        update_condition_text_match_counts(old_root)
                        old_root.save()
                        if old_root != ct:
                            update_condition_text_match_counts(ct)
                            ct.save()
                else:
                    # nothing has changed, no need to update count
                    pass
            else:
                ConditionTextMatch.objects.create(
                    parent=mode_of_inheritance_level,
                    condition_text=ct,
                    gene_symbol=gene_symbol,
                    mode_of_inheritance=mode_of_inheritance,
                    classification=classification
                )

                if update_counts:
                    update_condition_text_match_counts(ct)
                    ct.save()


def update_condition_text_match_counts(ct: ConditionText):
    # let's us count by doing one select out of the databse, rather than continually
    # selecting parents from the DB
    by_id: Dict[int, ConditionTextMatch] = dict()
    classification_related: List[ConditionTextMatch] = list()
    for ctm in ConditionTextMatch.objects.filter(condition_text=ct):
        by_id[ctm.id] = ctm
        if ctm.classification:
            classification_related.append(ctm)

    def check_hierarchy(ctm: ConditionTextMatch) -> bool:
        if ctm.is_valid:
            return True
        elif not ctm.is_blank:
            return False
        else:
            if parent := by_id.get(ctm.parent_id):
                return check_hierarchy(parent)
            return False

    invalid_count = 0
    for ctm in classification_related:
        if not check_hierarchy(ctm):
            invalid_count += 1

    ct.classifications_count = len(classification_related)
    ct.classifications_count_outstanding = invalid_count


@receiver(classification_post_publish_signal, sender=Classification)
def published(sender,
              classification: Classification,
              previously_published: ClassificationModification,
              newly_published: ClassificationModification,
              previous_share_level: ShareLevel,
              user: User,
              **kwargs):
    """
    Keeps condition_text_match in sync with the classifications when evidence changes
    """
    ConditionTextMatch.sync_condition_text_classification(newly_published)


@receiver(flag_comment_action, sender=Flag)
def check_for_discordance(sender, flag_comment: FlagComment, old_resolution: FlagResolution, **kwargs):
    """
    Keeps condition_text_match in sync with the classifications when withdraws happen/finish
    """
    flag = flag_comment.flag
    if flag.flag_type == flag_types.classification_flag_types.classification_withdrawn:
        cl: Classification
        if cl := Classification.objects.filter(flag_collection=flag.collection.id).first():
            ConditionTextMatch.sync_condition_text_classification(cl.last_published_version)
