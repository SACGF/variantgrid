from functools import reduce, cached_property
from typing import Any, Optional

from auditlog.models import AuditlogHistoryField
from auditlog.registry import auditlog
from django.contrib.auth.models import User
from django.db.models import CASCADE, QuerySet, PROTECT
from django.db import models
from django.db.models.enums import TextChoices
from django.utils.safestring import mark_safe
from django_extensions.db.models import TimeStampedModel
from annotation.models import ClinVarRecord
from classification.enums import OverlapStatus, TestingContextBucket, SpecialEKeys
from classification.models import ClassificationGrouping, EvidenceKeyMap, ConditionResolved, ClassificationResultValue
from classification.models.overlaps_enums import OverlapType, OverlapContributionStatus, TriageStatus, \
    OverlapEntrySourceTextChoices, EffectiveDateType, OverlapContributionChangeSource
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
#     Cached status of a contribution_status to the overlap.
#     Could be a reference to a classification or clinvar record
#     Useful to cache as the overall status is cached - so it's good to cache the working out
#     """
#     source: OverlapEntrySource
#     scv: Optional[str]
#     lab_id: Optional[int]
#     classification_grouping_id: Optional[int]
#     value: Optional[str]
#     # annoying thing about contribution_status is it takes a little bit of context knowledge to work out
#     contribution_status: OverlapContributionStatus = None #json_enum_encoder_for_text_choices(OverlapContributionStatus),
#     testing_context_bucket: TestingContextBucket = None #json_enum_encoder_for_text_choices(TestingContextBucket),
#     tumor_type_category: Optional[str] = None
#     # TODO do we want to keep date type somewhere?
#     effective_date: Optional[str] = None  # date str
#
#     @property
#     def _sort_obj(self):
#         # TODO do we need a sort order as part of value
#         return {self.source, self.testing_context_bucket, self.contribution_status, self.value, self.lab_id or 0}
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
    history = AuditlogHistoryField()

    source = models.TextField(choices=OverlapEntrySourceTextChoices.choices)
    scv = models.TextField(null=True, blank=True)  # could SCV change?
    allele = models.ForeignKey(Allele, null=True, blank=True, on_delete=CASCADE)
    classification_grouping = models.ForeignKey(ClassificationGrouping, null=True, blank=True, on_delete=CASCADE)
    value_type = models.TextField(choices=ClassificationResultValue.choices)
    value = models.TextField(null=True, blank=True)
    # annoying thing about contribution_status is it takes a little bit of context knowledge to work out
    contribution_status = models.TextField(choices=OverlapContributionStatus.choices)  # TODO rename to contribution_status
    testing_context_bucket = models.TextField(choices=TestingContextBucket.choices)
    tumor_type_category = models.TextField(null=True, blank=True)
    # TODO do we want to keep date type somewhere?
    effective_date = models.DateField(null=True, blank=True)
    effective_date_type = models.TextField(choices=EffectiveDateType.choices, default=EffectiveDateType.UNKNOWN)

    new_value = models.TextField(null=True, blank=True)  # TODO rename to amended value
    triage_status = models.TextField(max_length=1, choices=TriageStatus.choices, default=TriageStatus.PENDING)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._change_comment = None

    def set_change_comment(self, change_comment: str):
        self._change_comment = change_comment

    def get_additional_data(self) -> dict:
        return {"comment": self._change_comment}

    class Meta:
        unique_together = ('classification_grouping', 'value_type')

    @property
    def triage_status_obj(self) -> TriageStatus:
        return TriageStatus(self.triage_status)

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
        return OverlapContribution.pretty_value_for(self.value, self.value_type)

    @staticmethod
    def pretty_value_for(value: Optional[str], value_type: ClassificationResultValue):
        if not value:
            return "no-value"
        if value_type == ClassificationResultValue.ONC_PATH:
            return EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).pretty_value(value)
        elif value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            return EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(value)
        else:
            raise ValueError(f"Unsupported ValueType {value_type}")


auditlog.register(OverlapContribution)


class OverlapContributionLog(TimeStampedModel):
    overlap_contribution = models.ForeignKey(OverlapContribution, on_delete=CASCADE)
    change_source = models.TextField(choices=OverlapContributionChangeSource.choices)
    new_value = models.TextField(null=True, blank=True)  # TODO rename to amended value
    triage_status = models.TextField(max_length=1, choices=TriageStatus.choices, null=True, blank=True)
    comment = models.TextField(null=True, blank=True)
    user = models.ForeignKey(User, null=True, blank=True, on_delete=PROTECT)  # only used for triages

    # values changed by upload
    value = models.TextField(null=True, blank=True)
    effective_date = models.DateField(null=True, blank=True)
    effective_date_type = models.TextField(choices=EffectiveDateType.choices, default=EffectiveDateType.UNKNOWN)
    contribution_status = models.TextField(choices=OverlapContributionStatus.choices)

    @property
    def pretty_value(self) -> str:
        # FIXME remove duplication
        if not self.value:
            return "no-value"
        if self.overlap_contribution.value_type == ClassificationResultValue.ONC_PATH:
            return EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).pretty_value(self.value)
        elif self.overlap_contribution.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            return EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(self.value)
        else:
            raise ValueError(f"Unsupported ValueType {self.overlap_contribution.value_type}")

    @staticmethod
    def from_current_overlap_state(
            overlap_contribution: OverlapContribution,
            change_source: OverlapContributionChangeSource,
            comment: Optional[str] = None,
            user: Optional[User] = None):
        return OverlapContributionLog(
            overlap_contribution=overlap_contribution,
            change_source=change_source,
            new_value=overlap_contribution.new_value,
            triage_status=overlap_contribution.triage_status,
            comment=comment,
            user=user,
            value=overlap_contribution.value,
            effective_date=overlap_contribution.effective_date,
            effective_date_type=overlap_contribution.effective_date_type,
            contribution_status=overlap_contribution.contribution_status
        )


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
            for contribution in self.contributions.filter(contribution_status=OverlapContributionStatus.CONTRIBUTING).all():
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
            if entry.contribution_status == OverlapContributionStatus.CONTRIBUTING:
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

    @property
    def user_should_action(self) -> bool:
        match self:
            case TriageNextStep.AWAITING_YOUR_TRIAGE: return True
            case TriageNextStep.AWAITING_YOUR_AMEND: return True
            case TriageNextStep.TO_DISCUSS: return True
            case _: return False
    
    @property
    def icon(self):
        match self:
            case TriageNextStep.AWAITING_YOUR_TRIAGE:
                return mark_safe('<i class="fa-solid fa-clock mr-1" style="opacity:0.6"></i>')
            case TriageNextStep.AWAITING_YOUR_AMEND:
                return mark_safe('<i class="fa-solid fa-square-pen mr-1" style="opacity:0.6"></i>')
            case TriageNextStep.TO_DISCUSS:
                return mark_safe('<i class="fa-solid fa-comments mr-1" style="opacity:0.6"></i>')
            case _: return ""


# this should be the model that links Contributions to Overlaps to reduce redundancy
class OverlapContributionSkew(TimeStampedModel):
    overlap = models.ForeignKey(Overlap, on_delete=CASCADE)
    contribution = models.ForeignKey(OverlapContribution, on_delete=CASCADE)
    skew_perspective = models.TextField(max_length=1, choices=TriageNextStep.choices, default=TriageNextStep.PENDING_CALCULATION)
