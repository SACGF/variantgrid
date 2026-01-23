from functools import reduce, cached_property
from typing import Any, Optional

from django.contrib.auth.models import User
from django.db.models import CASCADE, QuerySet
from django.db import models
from django.db.models.enums import TextChoices
from django_extensions.db.models import TimeStampedModel
from annotation.models import ClinVarRecord
from classification.enums import OverlapStatus, TestingContextBucket, SpecialEKeys
from classification.models import ClassificationGrouping, EvidenceKeyMap, ConditionResolved, ClassificationResultValue
from classification.models.overlaps_enums import OverlapType, OverlapContributionStatus, TriageStatus, \
    OverlapEntrySourceTextChoices
from library.utils import first
from ontology.models import OntologyTerm
from snpdb.models import Allele, Lab


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
    scv = models.TextField(null=True, blank=True)  # could SCV change?
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

    @property
    def overlaps(self) -> QuerySet['Overlap']:
        return Overlap.objects.filter(pk__in=self.overlapcontributionskew_set.values_list('overlap_id'))

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

    @property
    def triage_status(self) -> TriageStatus:
        # TODO if we make triaging more complicated than just a date
        # we might want to
        if classification_grouping := self.classification_grouping:
            triage = classification_grouping.triage_for(self.value_type)
            return triage.triage_status
        else:
            return TriageStatus.NON_INTERACTIVE_THIRD_PARTY

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
    testing_context = models.TextField(max_length=1, choices=TestingContextBucket.choices, null=True, blank=True)
    tumor_type_category = models.TextField(null=True, blank=True)  # condition isn't always relevant
    overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_CONTRIBUTIONS.value)
    valid = models.BooleanField(default=False)  # if it's cross context but only has contributions from 1 context, or if it's NO_SUBMITTERS it shouldn't be valid

    # have to cache the values
    # contributions = models.ManyToManyField(OverlapContribution)
    @property
    def contributions(self) -> QuerySet[OverlapContribution]:
        return OverlapContribution.objects.filter(pk__in=self.overlapcontributionskew_set.values_list('contribution', flat=True))

    class Meta:
        indexes = [models.Index(fields=['overlap_type']), models.Index(fields=['value_type']), models.Index(fields=['allele'])]
        # TODO, we could put lab back in for ClinVar type so we can have this unique
        # unique_together = ('overlap_type', 'allele', 'value_type', 'testing_contexts', 'tumor_type_category', 'lab')

    @property
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
    def overlap_status_label(self):
        label = self.overlap_status_obj.label
        if self.overlap_type == OverlapType.CROSS_CONTEXT:
            match self.overlap_status_obj:
                case OverlapStatus.MAJOR_DIFFERENCES: return "Difference"
                case OverlapStatus.MEDICALLY_SIGNIFICANT: return "Medically significant difference"
                case _: return self.overlap_status_obj.label
        return self.overlap_status_obj.label

    @property
    def testing_contexts_objs(self) -> list[TestingContextBucket]:
        if testing_context := self.testing_context:
            return [TestingContextBucket(testing_context)]
        else:
            testing_contexts = set()
            for contribution in self.contributions.filter(contribution=OverlapContributionStatus.CONTRIBUTING).all():
                testing_contexts.add(TestingContextBucket(contribution.testing_context_bucket))

            return list(sorted(testing_contexts))

    @property
    def testing_context_obj(self) -> TestingContextBucket:
        if testing_contexts_objs := self.testing_contexts_objs:
            if len(testing_contexts_objs):
                return first(testing_contexts_objs)
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
            return list(EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).sort_values(relevant_values))[::-1]
        elif self.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            somatic_clin_sig_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
            return [somatic_clin_sig_e_key.pretty_value(val) for val in somatic_clin_sig_e_key.sort_values(relevant_values)][::-1]
        else:
            return []


class TriageNextStep(TextChoices):
    PENDING_CALCULATION = "X", "Pending Calculation"
    UNANIMOUSLY_COMPLEX = "C", "Unanimously Complex"
    AWAITING_YOUR_TRIAGE = "T", "Awaiting Your Triage"
    AWAITING_YOUR_AMEND = "A", "Awaiting Your Amendment"
    AWAITING_OTHER_LAB = "O", "Awaiting Other Lab"
    TO_DISCUSS = "D", "To Discuss"
    NOT_INVOLVED = "N", "Not Involved"


# this should be the model that links Contributions to Overlaps to reduce redundancy
class OverlapContributionSkew(TimeStampedModel):
    overlap = models.ForeignKey(Overlap, on_delete=CASCADE)
    contribution = models.ForeignKey(OverlapContribution, on_delete=CASCADE)
    skew_perspective = models.TextField(max_length=1, choices=TriageNextStep.choices, default=TriageNextStep.PENDING_CALCULATION)

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

    @cached_property
    def contribution(self) -> OverlapContribution:
        return OverlapContribution.objects.filter(classification_grouping=self.classification_grouping, value_type=self.result_value_type).get()

    @property
    def triage_status_obj(self):
        return TriageStatus(self.triage_status)

    class Meta:
        unique_together = ('classification_grouping', 'result_value_type')


class ClassificationGroupingValueTriageHistory(TimeStampedModel):
    triage = models.ForeignKey(ClassificationGroupingValueTriage, on_delete=models.CASCADE)
    new_value = models.TextField(null=True, blank=True)
    triage_status = models.TextField(max_length=1, choices=TriageStatus.choices, default=TriageStatus.PENDING)
    comment = models.TextField(null=True, blank=True)
    user = models.ForeignKey(User, null=False, blank=False, on_delete=models.PROTECT)
    state_data = models.JSONField(null=True, blank=True)  # For caching state of the Overlaps at the time

    @property
    def result_value_type(self):
        return self.triage.result_value_type

    @property
    def triage_status_obj(self):
        return TriageStatus(self.triage_status)