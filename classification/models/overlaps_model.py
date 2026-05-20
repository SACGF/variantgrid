from functools import reduce, cached_property
from typing import Any, Optional, Tuple
from auditlog.models import AuditlogHistoryField
from auditlog.registry import auditlog
from django.db.models import CASCADE, QuerySet, SET_NULL
from django.db import models
from django.db.models.enums import TextChoices, IntegerChoices
from django.urls import reverse
from django.utils.safestring import mark_safe
from django_extensions.db.models import TimeStampedModel
from annotation.models import ClinVarRecord, EffectiveDate
from classification.enums import OverlapStatus, TestingContextBucket, SpecialEKeys, TestingContextFull
from classification.models import ClassificationGrouping, EvidenceKeyMap, ConditionResolved, ClassificationResultValue
from classification.models.overlaps_enums import OverlapType, OverlapContributionStatus, OverlapEntrySourceTextChoices, \
    TriageState, TriageComment
from library.utils import first, AuditUtils
from library.utils.database_utils import TextFieldChoices, JSONDataclassField, IntegerFieldChoices
from ontology.models import OntologyTerm
from snpdb.models import Allele, Lab


class OverlapContribution(TimeStampedModel):
    history = AuditlogHistoryField()

    source = models.TextField(choices=OverlapEntrySourceTextChoices.choices)
    allele = models.ForeignKey(Allele, null=True, blank=True, on_delete=CASCADE)
    classification_grouping = models.ForeignKey(ClassificationGrouping, null=True, blank=True, on_delete=SET_NULL)
    scv = models.TextField(null=True, blank=True)  # could SCV change?

    value_type = models.TextField(choices=ClassificationResultValue.choices)
    value = models.TextField(null=True, blank=True)
    # annoying thing about contribution_status is it takes a little bit of context knowledge to work out

    contribution_status = TextFieldChoices(choices_type=OverlapContributionStatus)    # type: OverlapContributionStatus
    testing_context_bucket = models.TextField(choices=TestingContextBucket.choices)
    tumor_type_category = models.TextField(null=True, blank=True)

    effective_date = JSONDataclassField(dataclass_type=EffectiveDate, null=False, blank=False, default=EffectiveDate.default_json)  # type: EffectiveDate
    triage_state = JSONDataclassField(dataclass_type=TriageState, null=False, blank=False, default=TriageState.default_json)  # type: TriageState
    comment = JSONDataclassField(dataclass_type=TriageComment, null=False, blank=False, default=TriageComment.default_json)  # type: TriageComment

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # self._change_comment = None

    @cached_property
    def last_comment(self):
        return AuditUtils.last_change_for(self, "comment", is_json=True, parser=lambda x: TriageComment.from_dict(x))

    @property
    def testing_context_full(self) -> TestingContextFull:
        return TestingContextFull(
            testing_context_bucket=self.testing_context_bucket_obj,
            tumor_type_category=self.tumor_type_category
        )

    @property
    def effective_value(self):
        return self.triage_state.amend_value or self.value

    class Meta:
        unique_together = ('classification_grouping', 'value_type')

    @property
    def overlaps(self) -> QuerySet['Overlap']:
        return Overlap.objects.filter(pk__in=self.overlapcontributionskew_set.values_list('overlap_id'))

    @property
    def label(self):
        return f"{self.allele} {self.testing_context_full} {self.value_type}"

    def __str__(self):
        return f"{self.pk} {self.source} {self.value}"

    @property
    def testing_context_bucket_obj(self) -> TestingContextBucket:
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
                    terms = []
                    plain_texts = []
                    for condition_str in condition_strs:
                        try:
                            term = OntologyTerm.get_or_stub(condition_str)
                            terms.append(term)
                        except ValueError:
                            plain_texts.append(condition_str)
                    return ConditionResolved.from_uncounted_terms(terms=terms, plain_text_terms=plain_texts)
        return None

    @property
    def pretty_value(self) -> str:
        return OverlapContribution.pretty_value_for(self.value, self.value_type)

    @property
    def value_sort_index(self):
        if self.value_type == ClassificationResultValue.ONC_PATH:
            return EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).classification_sorter_value(self.value)
        elif self.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            return EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).classification_sorter_value(self.value)
        else:
            return 0

    def __lt__(self, other):
        if value_sort_diff := self.value_sort_index - other.value_sort_index:
            return value_sort_diff < 0
        if self.lab is None or other.lab is None:
            if self.lab is None and other.lab is None:
                return False
            if self.lab is None:
                return True
            else:
                return False
        return self.lab < other.lab

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


class Overlap(TimeStampedModel):
    """
    Overlap is made by composition as making a separate model for each overlap type added a lot of overhead for just isolating a few fields
    """
    history = AuditlogHistoryField()

    overlap_type = models.TextField(choices=OverlapType.choices)
    value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
    allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)  # might be blank for gene symbol wide
    testing_context_bucket = models.TextField(max_length=1, choices=TestingContextBucket.choices, null=True, blank=True)
    tumor_type_category = models.TextField(null=True, blank=True)  # condition isn't always relevant
    overlap_status = IntegerFieldChoices(choices_type=OverlapStatus, default=OverlapStatus.NO_CONTRIBUTIONS.value)  # type:OverlapStatus
    overlap_status_change_timestamp = models.DateTimeField(null=True, blank=True)
    valid = models.BooleanField(default=False)  # if it's cross context but only has contributions from 1 context, or if it's NO_SUBMITTERS it shouldn't be valid

    def get_absolute_url(self):
        return reverse('overlap_3', kwargs={"overlap_id": self.pk})

    # have to cache the values
    # contributions = models.ManyToManyField(OverlapContribution)
    @property
    def contributions(self) -> QuerySet[OverlapContribution]:
        return OverlapContribution.objects.filter(
            contribution_status=OverlapContributionStatus.CONTRIBUTING,
            pk__in=self.overlapcontributionskew_set.values_list('contribution', flat=True)
        ).select_related("classification_grouping__lab__organization")

    @property
    def contributions_all(self) -> QuerySet[OverlapContribution]:
        # unlike contributions this will also return OverlapContributions that aren't currently contribution
        # as they may have contributed in the past
        return  OverlapContribution.objects.filter(
            pk__in=self.overlapcontributionskew_set.values_list('contribution', flat=True)
        ).select_related("classification_grouping__lab__organization")

    @cached_property
    def contributions_list(self) -> list[OverlapContribution]:
        return list(self.contributions.all())

    class Meta:
        indexes = [models.Index(fields=['overlap_type']), models.Index(fields=['value_type']), models.Index(fields=['allele'])]
        # TODO, we could put lab back in for ClinVar type so we can have this unique
        # unique_together = ('overlap_type', 'allele', 'value_type', 'testing_contexts', 'tumor_type_category', 'lab')

    @property
    def scope_description(self):
        # at what scope is this overlap for, in future there could be gene symbol wide scopes
        if allele := self.allele:
            return str(allele)
        else:
            return "Unknown"

    @property
    def value_type_label(self):
        value_type = ClassificationResultValue(self.value_type)
        if value_type == ClassificationResultValue.ONC_PATH:
            if self.testing_contexts_objs == [TestingContextBucket.GERMLINE]:
                return "Pathogenicity"
        return value_type.label

    @property
    def overlap_status_label(self):
        if self.overlap_type == OverlapType.CROSS_CONTEXT:
            match self.overlap_status:
                case OverlapStatus.MAJOR_DIFFERENCES: return "Difference"
                case OverlapStatus.MEDICALLY_SIGNIFICANT: return "Medically significant difference"
                case _: return self.overlap_status.label
        return self.overlap_status.label

    @property
    def testing_contexts_objs(self) -> list[TestingContextBucket]:
        if testing_context := self.testing_context_bucket:
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
    def testing_context_full(self) -> TestingContextFull:
        # TOD
        return TestingContextFull(testing_context_bucket=self.testing_context_bucket, tumor_type_category=self.tumor_type_category)

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
            if overlap_type != OverlapType.SINGLE_CONTEXT:
                parts.append(OverlapType(overlap_type).label)
        if value_type := self.value_type:
            parts.append(ClassificationResultValue(value_type).label)
        parts.append("-".join(t.label for t in self.testing_contexts_objs))
        if tumor_type_category := self.tumor_type_category:
            parts.append(tumor_type_category)
        return " ".join(parts) + f" : {OverlapStatus(self.overlap_status).label}"

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


class TriageNextStep(IntegerChoices):
    NOT_INVOLVED = 0, "Not Involved"
    PENDING_CALCULATION = 1, "Pending Calculation"
    AWAITING_OTHER_LAB = 2, "Awaiting Other Lab"
    AWAITING_YOUR_TRIAGE = 3, "Awaiting Your Triage"
    AWAITING_YOUR_TRIAGE_OTHERS_TRIAGED = 4, "Awaiting your Triage - others have triaged"
    AWAITING_YOUR_AMEND = 5, "Pending Your Amendment"
    UNANIMOUSLY_COMPLEX = 6, "Unanimously Complex"
    TO_DISCUSS = 7, "To Discuss"

    @property
    def user_should_action(self) -> bool:
        match self:
            case TriageNextStep.AWAITING_YOUR_TRIAGE: return True
            case TriageNextStep.AWAITING_YOUR_TRIAGE_OTHERS_TRIAGED: return True
            case TriageNextStep.AWAITING_YOUR_AMEND: return True
            case TriageNextStep.TO_DISCUSS: return True
            case _: return False
    
    @property
    def icon(self):
        match self:
            case TriageNextStep.AWAITING_YOUR_TRIAGE:
                return mark_safe('<i class="fa-solid fa-clock mr-1" style="opacity:0.6"></i>')
            case TriageNextStep.AWAITING_YOUR_TRIAGE_OTHERS_TRIAGED:
                #TODO show a more impatient clock
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
    # skew_perspective = models.TextField(max_length=1, choices=TriageNextStep.choices, default=TriageNextStep.PENDING_CALCULATION)
    next_step = IntegerFieldChoices(choices_type=TriageNextStep, default=TriageNextStep.PENDING_CALCULATION)

    def __str__(self):
        return f"overlap = {self.overlap}, contribution = {self.contribution}, perspective = {self.next_step}"


class OverlapDiscordanceNotification(TimeStampedModel):
    overlap = models.ForeignKey('Overlap', on_delete=CASCADE)
    old_status = IntegerFieldChoices(OverlapStatus)  # type:OverlapStatus
    new_status = IntegerFieldChoices(OverlapStatus)  # type:OverlapStatus
    notification_sent_date = models.DateTimeField(null=True, blank=True)

    @property
    def is_still_relevant(self):
        if self.old_status.is_discordant ^ self.new_status.is_discordant:
            return True
        # TODO is going from somewhat discordant to more discordant notification worthy?
        return False

    def __lt__(self, other):
        return self.overlap_id < other.overlap_id