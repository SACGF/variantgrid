from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Dict, Any
import re

from django.contrib.auth.models import User
from django.contrib.postgres.fields import ArrayField
from django.core.exceptions import ValidationError
from django.db.models import CASCADE, QuerySet, SET_NULL
from guardian.shortcuts import get_objects_for_user, assign_perm
from lazy import lazy
from model_utils.models import TimeStampedModel
from django.db import models

from classification.enums import SpecialEKeys
from classification.models import Classification, ClassificationModification
from classification.regexes import db_ref_regexes, DbRefRegex, DbRegexes
from genes.models import GeneSymbol
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import all_users_group
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

    def update_counts(self):
        global_okay = False
        by_gene: Dict[str, bool] = dict()
        by_classification: Dict[str, bool] = dict()

        if root := self.conditiontextmatch_set.filter(gene_symbol__isnull=True).first():
            global_okay = root.is_valid

        ctms_gene_symbols = self.conditiontextmatch_set.filter(gene_symbol__isnull=False, classification__isnull=True)
        for ctm in ctms_gene_symbols:
            by_gene[ctm.gene_symbol.symbol] = ctm.is_valid or (ctm.is_blank and global_okay)

        ctms_classifications = self.conditiontextmatch_set.filter(classification__isnull=False)
        for ctm in ctms_classifications:
            by_classification[ctm.gene_symbol] = ctm.is_valid or (ctm.is_blank and by_gene.get(ctm.gene_symbol.symbol))

        self.classifications_count = len(by_classification)
        self.classifications_count_outstanding = len([status for status in by_classification.values() if status is False])

    @staticmethod
    def normalize(text: str):
        text = text.lower()
        text = re.sub("[,;.]", " ", text)  # replace , ; . with spaces
        text = re.sub("[ ]{2,}", " ", text)  # replace multiple spaces with
        return text

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


class ConditionTextMatch(TimeStampedModel, GuardianPermissionsMixin):
    condition_text = models.ForeignKey(ConditionText, on_delete=CASCADE)

    gene_symbol = models.ForeignKey(GeneSymbol, on_delete=CASCADE, null=True, blank=True)
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
        if self.classification:
            h_list.append(self.classification)
        return h_list

    @property
    def leaf(self) -> Optional[Any]:
        if self.classification:
            return self.classification
        elif self.gene_symbol:
            return self.gene_symbol
        else:
            return None

    # resolved to this,
    condition_xrefs = ArrayField(models.TextField(blank=False), default=list)
    condition_multi_operation = models.CharField(max_length=1, choices=MultiCondition.choices, default=MultiCondition.NOT_DECIDED)

    @property
    def condition_matching_str(self) -> str:
        result = ''
        if self.condition_xrefs:
            result = ", ".join(self.condition_xrefs)
        if len(self.condition_xrefs) >= 2:
            if self.condition_multi_operation == MultiCondition.NOT_DECIDED:
                result += "; uncertain/co-occuring"
            elif self.condition_multi_operation == MultiCondition.CO_OCCURRING:
                result += "; co-occuring"
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
        unique_together = ("condition_text", "gene_symbol", "classification")

    @property
    def is_valid(self):
        return len(self.condition_xrefs) == 1 or \
               (len(self.condition_xrefs) > 1 and self.condition_multi_operation != MultiCondition.NOT_DECIDED)

    @property
    def is_blank(self):
        return len(self.condition_xrefs) == 0

    @lazy
    def parent(self) -> Optional['ConditionTextMatch']:
        qs = ConditionTextMatch.objects.filter(condition_text=self.condition_text)
        if self.classification:
            qs = qs.filter(classification__isnull=True).filter(gene_symbol=self.gene_symbol)
        elif self.gene_symbol:
            qs = qs.filter(gene_symbol__isnull=True)
        else:
            return None
        return qs.first()

    @lazy
    def children(self) -> QuerySet:
        qs = ConditionTextMatch.objects.filter(condition_text=self.condition_text)
        if self.classification:
            return ConditionTextMatch.objects.none()
        elif self.gene_symbol:
            return qs.filter(classification__isnull=False).filter(gene_symbol=self.gene_symbol).order_by('classification_id')
        else:
            return qs.filter(classification__isnull=True, gene_symbol__isnull=False).order_by('gene_symbol__symbol')

    @staticmethod
    def sync_all():
        cms = ClassificationModification.objects.filter(is_last_published=True, classification__withdrawn=False).select_related("classification", "classification__lab")
        cm: ClassificationModification
        for cm in cms:
            ConditionTextMatch.sync_condition_text_classification(cm=cm)

    @staticmethod
    def sync_condition_text_classification(cm: ClassificationModification):
        classification = cm.classification

        if classification.withdrawn:
            ConditionTextMatch.objects.filter(classification=classification).delete()
            return

        lab = classification.lab
        gene_str = cm.get(SpecialEKeys.GENE_SYMBOL)
        if gene_symbol := GeneSymbol.objects.filter(symbol=gene_str).first():
            raw_condition_text = cm.get(SpecialEKeys.CONDITION) or ""
            normalized = ConditionText.normalize(raw_condition_text)
            ct: ConditionText
            ct, ct_is_new = ConditionText.objects.get_or_create(normalized_text=normalized, lab=lab)

            # if condition text has changed, remove the old entries
            ConditionTextMatch.objects.filter(classification=classification).exclude(condition_text=ct).delete()

            matches = db_ref_regexes.search(text=normalized)
            ConditionTextMatch.objects.get_or_create(
                condition_text=ct,
                gene_symbol=None,
                classification=None,
                defaults={"condition_xrefs": [match.id_fixed for match in matches]}
            )

            ConditionTextMatch.objects.get_or_create(
                condition_text=ct,
                gene_symbol=gene_symbol,
                classification=classification
            )
            ConditionTextMatch.objects.get_or_create(
                condition_text=ct,
                gene_symbol=gene_symbol,
                classification=None
            )

            ct.update_counts()
            ct.save()