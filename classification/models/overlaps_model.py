from dataclasses import dataclass
from enum import StrEnum
from functools import reduce, cached_property
from typing import Any, TypedDict, Optional

from dataclasses_json import dataclass_json
from django.contrib.postgres.fields import ArrayField
from django.db.models import TextChoices, JSONField, CASCADE
from django.db import models
from django_extensions.db.models import TimeStampedModel

from annotation.models import ClinVarRecord
from classification.enums import OverlapStatus, TestingContextBucket, SpecialEKeys, AlleleOriginBucket
from classification.models import ClassificationGrouping, EvidenceKeyMap, ConditionResolved
from library.utils import first
from ontology.models import OntologyTerm
from snpdb.models import Allele, Lab


class OverlapType(TextChoices):
    SINGLE_CONTEXT = "context", "Single Context"
    CROSS_CONTEXT = "cross", "Cross Context"
    CLINVAR_EXPERT_PANEL = "clinvar", "ClinVar Expert Panel"
    # CROSS_GERMLINE_NON_CANCER = "cross_germ_non_cancer", "Cross Germline / Somatic Non-Cancer"
    # CROSS_GERMLINE_HAEM = "cross_germ_haem", "Cross Germline / Somatic Haem"
    # CROSS_GERMLINE_SOLID_TUMOR = "cross_germ_solid_tumor", "Cross Germline / Solid Tumour"
    # CROSS_HAEM_SOLID_TUMOR = "cross_haem_solid_tumor", "Cross Haem / Solid Tumour"

    @property
    def priority_order(self) -> int:
        match self:
            case OverlapType.SINGLE_CONTEXT: return 1
            case OverlapType.CLINVAR_EXPERT_PANEL: return 2
            case OverlapType.CROSS_CONTEXT: return 3

    def __lt__(self, other):
        return self.priority_order < other.priority_order


class ClassificationResultValue(TextChoices):
    # FIXME should be called value type
    ONC_PATH = "O", "Onco-Path"
    CLINICAL_SIGNIFICANCE = "S", "Clinical significance"

    @property
    def priority_order(self) -> int:
        match self:
            case ClassificationResultValue.ONC_PATH: return 1
            case ClassificationResultValue.CLINICAL_SIGNIFICANCE: return 2


class OverlapContributionStatus(TextChoices):
    PENDING_CALCULATION = "P", "Pending Calculation"
    CONTRIBUTING = "C", "Contributing"
    NOT_SHARED = "N", "Not-shared"
    NO_VALUE = "X", "No value"
    NON_COMPARABLE_VALUE = "Z", "Non-comparable value"  # e.g. Risk Factor


class TriageStatus(TextChoices):
    PENDING = "P", "Pending Triage"
    REVIEWED_WILL_FIX = "F", "Will Amend"
    REVIEWED_WILL_DISCUSS = "D", "For Joint Discussion"
    REVIEWED_SATISFACTORY = "R", "Confident in Classification"
    COMPLEX = "X", "Low Penetrance/Risk Allele etc"


class OverlapEntrySourceTextChoices(TextChoices):
    CLASSIFICATION = "CLASS", "CLASSIFICATION"
    CLINVAR = "CLIN", "CLINVAR"


# class OverlapEntrySource(StrEnum):
#     """
#     Is this entry referencing a classification within this VariantGrid system or
#     """
#     CLASSIFICATION = "CLASSIFICATION"
#     CLINVAR = "CLINVAR"
#
#
# @dataclass_json
# @dataclass
# class OverlapEntry:
#     """
#     Cached status of a contribution to the overlap.
#     Could be a reference to a classification or clinvar record
#     Useful to cache as the overall status is cached - so it's good to cache the working out
#     """
#     source: OverlapEntrySource
#     scv: Optional[str]
#     lab_id: Optional[int]
#     classification_grouping_id: Optional[int]
#     value: Optional[str]
#     # annoying thing about contribution is it takes a little bit of context knowledge to work out
#     contribution: OverlapContributionStatus = None #json_enum_encoder_for_text_choices(OverlapContributionStatus),
#     testing_context_bucket: TestingContextBucket = None #json_enum_encoder_for_text_choices(TestingContextBucket),
#     tumor_type_category: Optional[str] = None
#     # TODO do we want to keep date type somewhere?
#     effective_date: Optional[str] = None  # date str
#
#     @property
#     def _sort_obj(self):
#         # TODO do we need a sort order as part of value
#         return {self.source, self.testing_context_bucket, self.contribution, self.value, self.lab_id or 0}
#
#     def __lt__(self, other):
#         return self._sort_obj < other._sort_obj
#
#     @property
#     def pretty_value(self) -> str:
#         if self.value and self.value.startswith("tier"):
#             return EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(self.value)
#         else:
#             return EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).pretty_value(self.value)
#
#     @cached_property
#     def lab(self) -> Optional[Lab]:
#         if lab_id := self.lab_id:
#             return Lab.objects.get(pk=lab_id)
#         return None


class OverlapContribution(TimeStampedModel):
    source = models.TextField(choices=OverlapEntrySourceTextChoices.choices)
    scv = models.TextField(null=True, blank=True) # could SCV change?
    allele = models.ForeignKey(Allele, null=True, blank=True, on_delete=CASCADE)
    classification_grouping = models.ForeignKey(ClassificationGrouping, null=True, blank=True, on_delete=CASCADE)
    value_type = models.TextField(choices=ClassificationResultValue.choices)
    value = models.TextField(null=True, blank=True)
    # annoying thing about contribution is it takes a little bit of context knowledge to work out
    contribution = models.TextField(choices=OverlapContributionStatus.choices)
    testing_context_bucket = models.TextField(choices=TestingContextBucket.choices)
    tumor_type_category = models.TextField(null=True, blank=True)
    # TODO do we want to keep date type somewhere?
    effective_date = models.DateField(null=True, blank=True)

    def __str__(self):
        return f"{self.pk} {self.source} {self.value}"

    @property
    def testing_context_bucket_obj(self):
        return TestingContextBucket(self.testing_context_bucket)

    @property
    def lab(self) -> Optional[Lab]:
        if classification_grouping := self.classification_grouping:
            return classification_grouping.lab
        return None

    @cached_property
    def conditions(self) -> Optional[ConditionResolved]:
        # TODO, should this be cached?
        if classification_grouping := self.classification_grouping:
            return ConditionResolved.from_dict(classification_grouping.conditions)
        elif scv := self.scv:
            if record := ClinVarRecord.objects.filter(record_id=scv).first():
                if condition_strs := record.conditions:
                    return ConditionResolved.from_uncounted_terms(terms=[OntologyTerm.get_or_stub(condition_str) for condition_str in condition_strs])
        return None

    @property
    def pretty_value(self) -> str:
        if not self.value:
            return "no-value"
        if self.value_type == ClassificationResultValue.ONC_PATH:
            return EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).pretty_value(self.value)
        elif self.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            return EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(self.value)
        else:
            raise ValueError(f"Unsupported ValueType {self.value_type}")


class Overlap(TimeStampedModel):
    """
    Overlap is made by composition as making a separate model for each overlap type added a lot of overhead for just isolating a few fields
    """
    overlap_type = models.TextField(choices=OverlapType.choices)
    value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
    allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)  # might be blank for gene symbol wide
    testing_contexts = ArrayField(models.TextField(max_length=1, choices=TestingContextBucket.choices), null=True, blank=True)
    tumor_type_category = models.TextField(null=True, blank=True)  # condition isn't always relevant
    overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_CONTRIBUTIONS.value)
    valid = models.BooleanField(default=False)  # if it's cross context but only has contributions from 1 context, or if it's NO_SUBMITTERS it shouldn't be valid

    # have to cache the values
    contributions = models.ManyToManyField(OverlapContribution)


    class Meta:
        indexes = [models.Index(fields=['overlap_type']), models.Index(fields=['value_type']), models.Index(fields=['allele'])]
        # TODO, we could put lab back in for ClinVar type so we can have this unique
        # unique_together = ('overlap_type', 'allele', 'value_type', 'testing_contexts', 'tumor_type_category', 'lab')

    def overlap_status_obj(self) -> OverlapStatus:
        return OverlapStatus(self.overlap_status)

    @property
    def value_type_label(self):
        value_type = ClassificationResultValue(self.value_type)
        if value_type == ClassificationResultValue.ONC_PATH:
            if self.testing_contexts_objs == [TestingContextBucket.GERMLINE]:
                return "Pathogenicity"
        return value_type.label

    @property
    def testing_contexts_objs(self) -> list[TestingContextBucket]:
        if testing_contexts := self.testing_contexts:
            return list(sorted([TestingContextBucket(t) for t in testing_contexts]))
        else:
            return []

    @property
    def testing_context(self) -> TestingContextBucket:
        if len(self.testing_contexts_objs) == 1:
            return first(self.testing_contexts_objs)
        else:
            raise ValueError("Overlap has multiple testing contexts")


    @property
    def priority_order(self) -> Any:
        return (
            self.allele_id,
            OverlapType(self.overlap_type).priority_order,
            reduce(lambda x, y: x*100 + y.priority_order, self.testing_contexts_objs, 0),
            ClassificationResultValue(self.value_type).priority_order,
            self.tumor_type_category or "",
            # self.lab_id if self.lab_id else 0,
        )

    def __lt__(self, other):
        return self.priority_order < other.priority_order

    def __str__(self):
        parts = []
        if not self.valid:
            parts.append("NOT IMPORTANT OVERLAP:")
        if allele := self.allele:
            parts.append(f"{allele:CA}")
        # if lab := self.lab:
        #     parts.append(str(lab))
        if overlap_type := self.overlap_type:
            parts.append(OverlapType(overlap_type).name)
        if value_type := self.value_type:
            parts.append(ClassificationResultValue(value_type).name)
        parts.append("-".join(t.name for t in self.testing_contexts_objs))
        if tumor_type_category := self.tumor_type_category:
            parts.append(tumor_type_category)
        return " ".join(parts) + f" : {OverlapStatus(self.overlap_status).name}"

    def relevant_values(self) -> list[str]:
        relevant_values = set()
        for entry in self.contributions.all():
            if entry.contribution == OverlapContributionStatus.CONTRIBUTING:
                relevant_values.add(entry.value)
        # for cg in self.classificationgroupingoverlapcontribution_set.filter(contribution_status=OverlapContributionStatus.CONTRIBUTING):
        #     # FIXME check triages
        #     # FIXME should this entire method
        #     if self.value_type == ClassificationResultValue.ONC_PATH:
        #         if value := cg.classification_grouping.latest_cached_summary.get("pathogenicity", {}).get("classification"):
        #             relevant_values.add(value)
        #     elif self.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
        #         if value := cg.classification_grouping.latest_cached_summary.get("somatic", {}).get("clinical_significance"):
        #             relevant_values.add(value)

        if self.value_type == ClassificationResultValue.ONC_PATH:
            return list(EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).sort_values(relevant_values))
        elif self.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            somatic_clin_sig_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
            return [somatic_clin_sig_e_key.pretty_value(val) for val in somatic_clin_sig_e_key.sort_values(relevant_values)]
        else:
            return []


#
# @dataclass
# class OverlapSkew:
#     classification_grouping_id: int
#     overlap: Overlap
#     other_testing_contexts: list[TestingContextBucket]
#     other_context_entries: list[OverlapEntry]
#
#     @property
#     def overlap_entries(self) -> list[OverlapEntry]:
#         return self.other_context_entries
#
#     @staticmethod
#     def skew_versus(overlap: Overlap, classification_grouping_id: int) -> 'OverlapSkew':
#         matching_entry: Optional[OverlapEntry] = None
#         other_entries: list[OverlapEntry] = []
#         for entry in overlap.overlap_entries:
#             if entry.classification_grouping_id == classification_grouping_id:
#                 matching_entry = entry
#             else:
#                 other_entries.append(entry)
#
#         if not matching_entry:
#             raise ValueError(f"{classification_grouping_id} is not part of this entry")
#
#         match overlap.overlap_type:
#             case OverlapType.SINGLE_CONTEXT | OverlapType.CLINVAR_EXPERT_PANEL:
#                 return OverlapSkew(
#                     classification_grouping_id,
#                     overlap,
#                     other_testing_contexts=[],
#                     other_context_entries=other_entries)
#             case OverlapType.CROSS_CONTEXT:
#                 other_contexts = set(overlap.testing_contexts_objs)
#                 other_contexts.remove(matching_entry.testing_context_bucket)
#                 return OverlapSkew(
#                     classification_grouping_id,
#                     overlap,
#                     other_testing_contexts=list(other_contexts)[0],
#                     other_context_entries=[entry for entry in overlap.overlap_entries if entry.testing_context_bucket in other_contexts]
#                 )
#             case _:
#                 raise ValueError(f"Unexpected overlap type {overlap.overlap_type}")


# class ClassificationGroupingOverlapContribution(TimeStampedModel):
#     classification_grouping = models.ForeignKey(ClassificationGrouping, on_delete=models.CASCADE)
#     overlap = models.ForeignKey(Overlap, on_delete=models.CASCADE)
#     contribution_status = models.TextField(max_length=1, choices=OverlapContributionStatus.choices, default=OverlapContributionStatus.PENDING_CALCULATION)
#
#     class Meta:
#         unique_together = ('classification_grouping', 'overlap')
#

class ClassificationGroupingValueTriage(TimeStampedModel):
    classification_grouping = models.ForeignKey(ClassificationGrouping, on_delete=models.CASCADE)
    result_value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
    new_value = models.TextField(null=True, blank=True)
    triage_status = models.TextField(max_length=1, choices=TriageStatus.choices, default=TriageStatus.PENDING)

    class Meta:
        unique_together = ('classification_grouping', 'result_value_type')
