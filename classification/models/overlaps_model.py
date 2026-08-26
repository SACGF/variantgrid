from functools import reduce, cached_property
from typing import Any, Optional
from auditlog.models import AuditlogHistoryField
from auditlog.registry import auditlog
from django.conf import settings
from django.db.models import CASCADE, QuerySet, SET_NULL, JSONField
from django.db import models
from django.db.models.enums import IntegerChoices
from django.dispatch.dispatcher import receiver
from django.urls import reverse
from django.utils.safestring import mark_safe
from django_extensions.db.models import TimeStampedModel
from annotation.models import ClinVarRecord
from annotation.models.data_enums import EffectiveDate
from classification.enums import OverlapStatus, TestingContextBucket, SpecialEKeys, TestingContextFull, TriageStatus, \
    OverlapOverrideStatus, OverlapState
from classification.models import ClassificationGrouping, EvidenceKeyMap, ConditionResolved, ClassificationResultValue
from classification.enums.overlaps_enums import OverlapType, OverlapContributionStatus, OverlapEntrySourceTextChoices, \
    TriageState, TriageComment
from genes.hgvs import HGVSComponents, HGVSDisplay
from library.guardian_utils import admin_bot
from library.preview_request import PreviewModelMixin, PreviewKeyValue
from library.utils import first, AuditUtils, AuditSingleChange
from library.utils.database_utils import TextFieldChoices, IntegerFieldChoices
from ontology.models import OntologyTerm
from review.models import ReviewableModelMixin, Review, review_detail_signal
from snpdb.models import Allele, Lab, GenomeBuild, LabLike, CLINVAR_EXPERT_PANEL_LAB

IN_REVIEW_VALUE = "in-review"


class OverlapContribution(TimeStampedModel):
    history = AuditlogHistoryField()

    source = models.TextField(choices=OverlapEntrySourceTextChoices.choices)
    allele = models.ForeignKey(Allele, null=True, blank=True, on_delete=CASCADE)
    classification_grouping = models.ForeignKey(ClassificationGrouping, null=True, blank=True, on_delete=SET_NULL)
    scv = models.TextField(null=True, blank=True)  # could SCV change?

    value_type = models.TextField(choices=ClassificationResultValue.choices)
    value = models.TextField(null=True, blank=True)

    contribution_status = TextFieldChoices(choices_type=OverlapContributionStatus)    # type: OverlapContributionStatus
    testing_context_bucket = models.TextField(choices=TestingContextBucket.choices)
    tumor_type_category = models.TextField(null=True, blank=True)

    effective_date = JSONField(null=False, blank=False, default=EffectiveDate.default_json)
    triage_state = JSONField(null=False, blank=False, default=TriageState.default_json)  # pending values are captured here
    comment = JSONField(null=False, blank=False, default=TriageComment.default_json)

    review_agreed_value = models.TextField(null=True, blank=True)
    """
    If this is not None, if it matches the effective_value, and every other OverlapContribution in the Overlap has their
    final reviewed value matching their 
    """

    def is_review_agreed_value_met(self) -> bool:
        """
        Was there a review where there was going to be continued discordance
        :return:
        """
        if review_agreed_value := self.review_agreed_value:
            return review_agreed_value == self.effective_value
        else:
            return False

    @property
    def lab_like(self) -> LabLike:
        if cg := self.classification_grouping:
            return cg.lab
        elif self.scv:
            return CLINVAR_EXPERT_PANEL_LAB
        else:
            raise ValueError("Cannot determine LabLike")

    @property
    def effective_date_obj(self) -> EffectiveDate:
        try:
            return EffectiveDate.from_dict(self.effective_date)
        except Exception as ex:
            return EffectiveDate(date=self.effective_date)

    @effective_date_obj.setter
    def effective_date_obj(self, value: EffectiveDate):
        self.effective_date = value.to_dict()

    @property
    def triage_state_obj(self) -> TriageState:
        return TriageState.from_dict(self.triage_state)

    @triage_state_obj.setter
    def triage_state_obj(self, value: TriageState):
        self.triage_state = value.to_dict()

    @property
    def comment_obj(self) -> TriageComment:
        return TriageComment.from_dict(self.comment)

    @comment_obj.setter
    def comment_obj(self, value: TriageComment):
        self.comment = value.to_dict()

    @cached_property
    def last_comment(self) -> AuditSingleChange[TriageComment]:
        if last_comment := AuditUtils.last_change_for(self, "comment", is_json=True, parser=lambda x: TriageComment.from_dict(x)):
            if last_comment.value.count != 0:
                return last_comment
        if not self.triage_state_obj.status.is_default:
            if last_triage_change := AuditUtils.last_change_for(self, "triage_state", parser=lambda x: None):
                return last_triage_change
        return None

    @property
    def testing_context_full(self) -> TestingContextFull:
        return TestingContextFull(
            testing_context_bucket=TestingContextBucket(self.testing_context_bucket_obj),
            tumor_type_category=self.tumor_type_category
        )

    @property
    def effective_value(self):
        return self.triage_state_obj.amend_value or self.value

    @property
    def pending_value(self) -> Optional[str]:
        if self.is_amending:
            return self.triage_state_obj.amend_value or IN_REVIEW_VALUE
        return None

    @property
    def is_amending(self):
        return self.triage_state_obj.status == TriageStatus.REVIEWED_WILL_FIX

    class Meta:
        unique_together = ('classification_grouping', 'value_type')

    @property
    def overlaps(self) -> QuerySet['Overlap']:
        return Overlap.objects.filter(pk__in=self.overlapcontributionskew_set.values_list('overlap_id'))

    @property
    def label(self):
        return f"{self.allele} {self.testing_context_full} {self.value_type}"

    def __str__(self):
        return f"{self.pk} {self.source} {self.lab_like} {self.value}"

    @property
    def testing_context_bucket_obj(self) -> TestingContextBucket:
        return TestingContextBucket(self.testing_context_bucket)

    @property
    def lab(self) -> Optional[Lab]:
        # DEPRECATED, use lab_like
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
    def pretty_effective_value(self):
        return OverlapContribution.pretty_value_for(self.effective_value, self.value_type)

    @property
    def value_sort_index(self):
        if self.value_type == ClassificationResultValue.ONC_PATH:
            return EvidenceKeyMap.cached_key(SpecialEKeys.ONC_PATH).classification_sorter_value(self.value)
        elif self.value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE:
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
        elif value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE:
            return EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE).pretty_value(value)
        else:
            raise ValueError(f"Unsupported ValueType {value_type}")


auditlog.register(OverlapContribution)


class Overlap(TimeStampedModel, ReviewableModelMixin, PreviewModelMixin):
    """
    Overlap is made by composition as making a separate model for each overlap type added a lot of overhead for just isolating a few fields
    """
    history = AuditlogHistoryField()

    overlap_type = models.TextField(choices=OverlapType.choices)
    value_type = TextFieldChoices(max_length=1, choices_type=ClassificationResultValue)  # type:ClassificationResultValue
    allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)  # might be blank for gene symbol wide
    testing_context_bucket = models.TextField(max_length=1, choices=TestingContextBucket.choices, null=True, blank=True)
    tumor_type_category = models.TextField(null=True, blank=True)  # condition isn't always relevant
    overlap_status = IntegerFieldChoices(choices_type=OverlapStatus, default=OverlapStatus.NO_CONTRIBUTIONS.value)  # type:OverlapStatus
    has_pending_values = models.BooleanField(default=False, blank=True)  # if

    overlap_override_status = IntegerFieldChoices(choices_type=OverlapOverrideStatus, default=OverlapOverrideStatus.NO_OVERRIDE)  # type:OverlapOverrideStatus

    overlap_max_ever_status = IntegerFieldChoices(choices_type=OverlapStatus, default=OverlapStatus.NO_CONTRIBUTIONS.value)  # type:OverlapStatus
    overlap_status_change_timestamp = models.DateTimeField(null=True, blank=True)
    valid = models.BooleanField(default=False)  # if it's cross context but only has contributions from 1 context, or if it's NO_SUBMITTERS it shouldn't be valid

    cached_overlap_state = models.JSONField(default=OverlapState.default_json)

    @property
    def cached_overlap_state_obj(self) -> OverlapState:
        print(self.cached_overlap_state)
        return OverlapState.from_dict(self.cached_overlap_state)

    @cached_overlap_state_obj.setter
    def cached_overlap_state_obj(self, value: OverlapState):
        self.cached_overlap_state = value.to_dict()

    @property
    def is_active_discordance(self):
        return self.overlap_status.is_discordant and not self.overlap_override_status

    @property
    def is_active_supported_discordance(self) -> bool:
        # TODO move OVERLAP_CLIN_SIG_ENABLED to a more core location
        from classification.services.overlap_calculator import OVERLAP_CLIN_SIG_ENABLED
        if self.overlap_type == OverlapType.CROSS_CONTEXT:
            return False
        if self.value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE and not OVERLAP_CLIN_SIG_ENABLED:
            return False
        return self.is_active_discordance

    def get_absolute_url(self):
        return reverse('overlap_3', kwargs={"overlap_id": self.pk})

    @property
    def derived_overlap_state(self):
        """
        Cached overlap_state is pure JSON, but dervied does also blend key values
        It's a bit redundant but cached is used to see if we need to send out notifications
        whereas the database values and derived are used for rendering
        """
        return OverlapState(
            status=self.overlap_status,
            has_pending_values=self.has_pending_values,
            override_status=self.overlap_override_status,
            lab_groups=list(sorted(cont.classification_grouping.lab.group_name for cont in self.contributions_list if cont.classification_grouping is not None))
        )

    @classmethod
    def preview_category(cls) -> str:
        return "Overlap"

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-record-vinyl"

    @classmethod
    def preview_enabled(cls) -> bool:
        return settings.DISCORDANCE_ENABLED

    @property
    def preview(self) -> 'PreviewData':
        summary_extra = []
        if not self.valid:
            summary_extra.append(PreviewKeyValue(value="Invalid Overlap - ignore", dedicated_row=True))
        else:
            summary_extra.append(PreviewKeyValue("Status", ""))
            if allele := self.allele:
                summary_extra.append(PreviewKeyValue("Allele", f"{allele:CA}"))
            if overlap_type := self.overlap_type:
                if overlap_type != OverlapType.SINGLE_CONTEXT:
                    summary_extra.append(PreviewKeyValue(value=OverlapType(overlap_type).label))
            if value_type := self.value_type:
                summary_extra.append(PreviewKeyValue(value=ClassificationResultValue(value_type).label))

        for c_hgvs in self.c_hgvses:
            summary_extra.append(
                PreviewKeyValue(key=f"{c_hgvs.genome_build} c.HGVS", value=str(c_hgvs), dedicated_row=True)
            )

        return self.preview_with(
            identifier=f"OV_{self.pk}",
            summary_extra=summary_extra
        )

    @cached_property
    def c_hgvses(self):
        c_hgvses = set()
        for entry in self.contributions_list:
            if cg := entry.classification_grouping:
                c_hgvses.add(cg.latest_allele_info.preferred_c_hgvs_obj())
        return list(sorted(c_hgvses))

    def c_hgvs(self, lab: Lab, genome_build: Optional[GenomeBuild] = None) -> HGVSDisplay:
        # if no genome_build provided, use the imported value
        # TODO, if there are multiple contributions from the same lab, should we get multiple c.HGVSs?
        lab_classification_grouping = self.contributions.filter(classification_grouping__lab=lab).first()
        if not lab_classification_grouping:
            # the lab doesn't actually have a horse in this game
            lab_classification_grouping = self.contributions.filter(classification_grouping__isnull=False).first()
        if not lab_classification_grouping:
            return HGVSDisplay(components=HGVSComponents())  # got nothing to work with in this overlap
        if genome_build:
            return lab_classification_grouping.classification_grouping.latest_allele_info.preferred_c_hgvs_obj(genome_build)
        else:
            return lab_classification_grouping.classification_grouping.latest_allele_info.imported_c_hgvs_obj

    def c_hgvs_all(self, genome_build: GenomeBuild) -> list[HGVSDisplay]:
        results = set()
        for contribution in self.contributions_list:
            if classification_grouping := contribution.classification_grouping:
                results.add(classification_grouping.latest_allele_info.preferred_c_hgvs_obj(genome_build))
        return list(sorted(results))

    # have to cache the values
    # contributions = models.ManyToManyField(OverlapContribution)
    @property
    def contributions(self) -> QuerySet[OverlapContribution]:
        return OverlapContribution.objects.filter(
            contribution_status=OverlapContributionStatus.CONTRIBUTING,
            pk__in=self.overlapcontributionskew_set.values_list('contribution', flat=True)
        ).select_related("classification_grouping__lab__organization")

    @property
    def reviewing_labs(self) -> set[Lab]:
        lab_ids = set(self.contributions.values_list('classification_grouping__lab', flat=True))
        return set(Lab.objects.filter(pk__in=lab_ids).all())

    def post_review_url(self, review: Review) -> str:
        return reverse('action_overlap_review', kwargs={"review_id": review.pk})

    @property
    def has_clinvar_expert_panel(self) -> bool:
        return self.contributions.filter(scv__isnull=False).exists()

    @property
    def contributions_all(self) -> QuerySet[OverlapContribution]:
        # unlike contributions this will also return OverlapContributions that aren't currently contribution
        # as they may have contributed in the past
        return OverlapContribution.objects.filter(
            pk__in=self.overlapcontributionskew_set.values_list('contribution', flat=True)
        ).select_related("classification_grouping__lab__organization")

    @cached_property
    def contributions_list(self) -> list[OverlapContribution]:
        return list(sorted(self.contributions.all()))

    @property
    def has_outstanding_triages(self) -> bool:
        return self.contributions.filter(triage_state__status=TriageStatus.PENDING).exists()

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
        value_type = self.value_type
        if value_type == ClassificationResultValue.ONC_PATH:
            if self.testing_contexts_objs == [TestingContextBucket.GERMLINE]:
                return "Pathogenicity"
            elif TestingContextBucket.GERMLINE not in self.testing_contexts_objs:
                return "Oncogenicity"
        else:
            return "Clinical Significance"
        return value_type.label

    @property
    def is_single_context(self):
        return self.overlap_type == OverlapType.SINGLE_CONTEXT

    @property
    def overlap_status_label(self):
        if self.overlap_type == OverlapType.CROSS_CONTEXT or self.value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE:
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
        return TestingContextFull(
            testing_context_bucket=TestingContextBucket(self.testing_context_bucket),
            tumor_type_category=self.tumor_type_category
        )

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
                relevant_values.add(entry.effective_value)
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
            return list(val.replace('_', '-') for val in EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE).sort_values(relevant_values))
        elif self.value_type == ClassificationResultValue.SOMATIC_CLINICAL_SIGNIFICANCE:
            tier_1_or_2 = False
            if "tier_1_or_2" in relevant_values:
                relevant_values.remove("tier_1_or_2")
                tier_1_or_2 = True
            somatic_clin_sig_e_key = EvidenceKeyMap.cached_key(SpecialEKeys.SOMATIC_CLINICAL_SIGNIFICANCE)
            sorted_values = somatic_clin_sig_e_key.sort_values(relevant_values)[::-1]
            if tier_1_or_2:
                sorted_values.append("tier_1_or_2")

            result = [somatic_clin_sig_e_key.pretty_value(val) for val in sorted_values]
            return result
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
                # TODO show a more impatient clock
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
    next_step = IntegerFieldChoices(choices_type=TriageNextStep, default=TriageNextStep.PENDING_CALCULATION)

    def __str__(self):
        return f"{self.contribution} ({self.next_step.name})"


class OverlapDiscordanceNotification(TimeStampedModel):
    overlap = models.ForeignKey('Overlap', on_delete=CASCADE)
    old_state = models.JSONField(default=OverlapState.default_json)
    new_state = models.JSONField(default=OverlapState.default_json)

    notification_sent_date = models.DateTimeField(null=True, blank=True)

    @property
    def old_state_obj(self):
        return OverlapState.from_dict(self.old_state)

    @old_state_obj.setter
    def old_state_obj(self, value: OverlapState):
        self.old_state = value.to_dict()

    @property
    def new_state_obj(self):
        return OverlapState.from_dict(self.new_state)

    @new_state_obj.setter
    def new_state_obj(self, value: OverlapState):
        self.new_state = value.to_dict()

    @property
    def is_still_relevant(self):
        return OverlapState.is_notify_relevant(self.old_state_obj, self.new_state_obj)

    def __lt__(self, other):
        return self.overlap_id < other.overlap_id
