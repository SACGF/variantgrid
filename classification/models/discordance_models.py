import typing
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Set, Optional, List, Dict

import django.dispatch
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models, transaction
from django.db.models.deletion import PROTECT, CASCADE
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from lazy import lazy

from classification.enums.classification_enums import SpecialEKeys, ClinicalSignificance
from classification.enums.discordance_enums import DiscordanceReportResolution, ContinuedDiscordanceReason
from classification.models.classification import ClassificationModification, Classification
from classification.models.clinical_context_models import ClinicalContext
from classification.models.flag_types import classification_flag_types
from flags.models.enums import FlagStatus
from flags.models.models import FlagComment
from genes.hgvs import CHGVS
from snpdb.genome_build_manager import GenomeBuildManager
from snpdb.models import Lab, GenomeBuild

discordance_change_signal = django.dispatch.Signal()  # args: "discordance_report"


class DiscordanceReport(TimeStampedModel):

    resolution = models.TextField(default=DiscordanceReportResolution.ONGOING, choices=DiscordanceReportResolution.CHOICES, max_length=1, null=True, blank=True)
    continued_discordance_reason = models.TextField(choices=ContinuedDiscordanceReason.CHOICES, max_length=1, null=True, blank=True)
    continued_discordance_text = models.TextField(null=True, blank=True)

    clinical_context = models.ForeignKey(ClinicalContext, on_delete=PROTECT)

    # report started date is a bit redundant compared to TimeStampModel's created
    # but as it has more of a business implication than a technical one I felt it
    # should be its own field
    report_closed_by = models.ForeignKey(User, on_delete=PROTECT, null=True)
    report_started_date = models.DateTimeField(auto_now_add=True)
    report_completed_date = models.DateTimeField(null=True, blank=True)

    cause_text = models.TextField(null=False, blank=True, default='')
    resolved_text = models.TextField(null=False, blank=True, default='')

    def get_absolute_url(self):
        return reverse('discordance_report', kwargs={"report_id": self.pk})

    @property
    def status(self) -> str:
        if self.resolution:
            return self.get_resolution_display()
        return "Ongoing"

    @property
    def is_important(self):
        return self.resolution is None or self.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE

    @property
    def is_active(self) -> bool:
        return self.resolution is None

    @property
    def resolution_text(self) -> str:
        if self.resolution == DiscordanceReportResolution.CONCORDANT:
            return 'Concordant'
        if self.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
            return 'Continued Discordance'
        if self.resolution == DiscordanceReportResolution.ONGOING:
            return 'In Discordance'
        return ''

    @transaction.atomic
    def close(self, expected_resolution: Optional[str] = None, cause_text: str = ''):
        """
        @param expected_resolution this will be calculated, but if provided a ValueError will be raised if it doesn't match
        @param cause_text reason this discordance is being closed
        """
        if self.clinical_context.is_discordant():
            self.resolution = DiscordanceReportResolution.CONTINUED_DISCORDANCE
        else:
            self.resolution = DiscordanceReportResolution.CONCORDANT

        if expected_resolution and expected_resolution != self.resolution:
            raise ValueError(f'Expected to close Discordance Report as {expected_resolution} but was {self.resolution}')

        self.resolved_text = cause_text
        self.report_completed_date = datetime.now()
        self.save()

        for drc in DiscordanceReportClassification.objects.filter(report=self):  # type: DiscordanceReportClassification
            drc.close()

        discordance_change_signal.send(DiscordanceReport, discordance_report=self)

    @transaction.atomic
    def unresolve_close(self, user: User, continued_discordance_reason: str, continued_discordance_text: str):
        self.report_closed_by = user
        self.continued_discordance_reason = continued_discordance_reason
        self.continued_discordance_text = continued_discordance_text
        self.close(expected_resolution=DiscordanceReportResolution.CONTINUED_DISCORDANCE, cause_text="Unable to resolve")

    @transaction.atomic
    def update(self, cause_text: str = ''):
        if not self.is_active:
            raise ValueError('Cannot update non active Discordance Report')

        if not self.id:
            self.save()

        existing_vms = set()
        existing_labs = set()
        for drc in DiscordanceReportClassification.objects.filter(report=self):
            existing_vms.add(drc.classification_original.classification_id)
            existing_labs.add(drc.classification_original.classification.lab)

        newly_added_labs: Set[Lab] = set()
        for vcm_id in self.clinical_context.classifications_qs.values_list('id', flat=True):
            if vcm_id in existing_vms:
                existing_vms.remove(vcm_id)
            else:
                vcm = ClassificationModification.objects.get(is_last_published=True, classification=vcm_id)
                if vcm.classification.lab not in existing_labs:
                    newly_added_labs.add(vcm.classification.lab)

                DiscordanceReportClassification(
                    report=self,
                    classification_original=vcm
                ).save()

        # contents of existing_vms are now presumably records that changed clinical groupings

        if not self.clinical_context.is_discordant():
            # hey we've reached concordance, let's close this
            # close will fire off a notification event
            self.close(expected_resolution=DiscordanceReportResolution.CONCORDANT, cause_text=cause_text)
        else:
            if newly_added_labs:  # change is significant
                discordance_change_signal.send(DiscordanceReport, discordance_report=self)

    def all_actively_involved_labs(self):
        all_lab_ids = set()
        for lab_id in DiscordanceReportClassification.objects.filter(report=self).values_list('classification_original__classification__lab', flat=True):
            all_lab_ids.add(lab_id)
        return list(Lab.objects.filter(id__in=all_lab_ids).all())

    @transaction.atomic
    def create_new_report(self, only_if_necessary: bool = True, cause: str = ''):
        if not self.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
            raise ValueError(f'Should only call create_new_report from the latest report, only if resolution = {DiscordanceReportResolution.CONTINUED_DISCORDANCE}')

        report = None
        if (not only_if_necessary) or self.has_significance_changed():
            report = DiscordanceReport(clinical_context=self.clinical_context, cause_text=cause)
            report.save()
            # if we're creating a new discordance report after a "Continued Discordance" we want to grab the state from as it was at the time
            # the previous report was closed, not the current date
            for last_id in DiscordanceReportClassification.objects.filter(
                    report=self,
                    clinical_context_final=self.clinical_context
                ).values_list('classification_final', flat=True):

                DiscordanceReportClassification(
                    report=self,
                    classification_original=ClassificationModification.objects.get(pk=last_id)
                ).save()

            # the report might even auto-close itself if the change brought it into concordance
            report.update(cause_text=cause)
        return report

    def has_significance_changed(self):
        existing_vms = {}
        for dr in DiscordanceReportClassification.objects.filter(report=self):
            if dr.clinical_context_final == self.clinical_context:
                existing_vms[dr.classification_final.classification.id] = dr.classification_final.get(SpecialEKeys.CLINICAL_SIGNIFICANCE, None)

        for vcm in self.clinical_context.classification_modifications:
            vc_id = vcm.classification.id
            if vc_id not in existing_vms:
                return True  # new classification, previous state is not relevant
            if existing_vms[vc_id] != vcm.get(SpecialEKeys.CLINICAL_SIGNIFICANCE):
                return True  # clinical significance has changed
            existing_vms.pop(vc_id, None)

        if existing_vms:
            return True
        return False

    @staticmethod
    def latest_report(clinical_context: ClinicalContext) -> 'DiscordanceReport':
        return DiscordanceReport.objects.filter(clinical_context=clinical_context).order_by('-created').first()

    @staticmethod
    def update_latest(clinical_context: ClinicalContext, cause: str):
        latest_report = DiscordanceReport.latest_report(clinical_context=clinical_context)
        if latest_report:
            if latest_report.is_active:
                latest_report.update(cause_text=cause)
                return latest_report
            if latest_report.resolution == DiscordanceReportResolution.CONCORDANT:
                # most recent report ended with discordance, can start fresh
                pass
            elif latest_report.resolution == DiscordanceReportResolution.CONTINUED_DISCORDANCE:
                # create a new report if data has changed significantly
                return latest_report.create_new_report(only_if_necessary=True, cause=f'Change after previous report marked as continued discordance - {cause}')

        if clinical_context.is_discordant():
            report = DiscordanceReport(clinical_context=clinical_context, cause_text=cause)
            report.update(cause_text=cause)
            return report

        # no report has been made
        return None

    def actively_discordant_classification_ids(self) -> Set[Classification]:
        if not self.is_active:
            return set()
        classifications = set()
        for dr in DiscordanceReportClassification.objects.filter(report=self):
            if dr.clinical_context_effective == self.clinical_context and not dr.withdrawn_effective:
                classifications.add(dr.classification_original.classification.id)
        return classifications

    @lazy
    def all_classification_modifications(self) -> List[ClassificationModification]:
        vcms: List[ClassificationModification] = list()
        for dr in DiscordanceReportClassification.objects.filter(report=self):
            vcms.append(dr.classification_effective)
        return vcms

    def all_c_hgvs(self, genome_build: Optional[GenomeBuild] = None) -> List[CHGVS]:
        if not genome_build:
            genome_build = GenomeBuildManager.get_current_genome_build()
        c_hgvs = set()
        for cm in self.all_classification_modifications:
            c_hgvs.add(cm.c_hgvs_best(genome_build))
        return sorted(c_hgvs)

    @dataclass
    class DiscordanceLabSummary:
        lab: Lab
        values: Set[str]

        def __lt__(self, other):
            return self.lab < other.lab

    def labs_to_classification(self) -> List[DiscordanceLabSummary]:
        lab_to_class: Dict[Lab, Set[str]] = defaultdict(set)
        for drc in DiscordanceReportClassification.objects.filter(report=self):
            lab_to_class[drc.classification_effective.classification.lab].add(drc.classification_effective.get(SpecialEKeys.CLINICAL_SIGNIFICANCE))

        from classification.models import EvidenceKeyMap
        clin_sig = EvidenceKeyMap.cached_key(SpecialEKeys.CLINICAL_SIGNIFICANCE)

        return sorted([DiscordanceReport.DiscordanceLabSummary(k, clin_sig.sort_values(v)) for k, v in lab_to_class.items()])


class DiscordanceAction:

    def __init__(self, code: int, short_description: str, long_description: str = None, no_longer_considered: bool = False):
        self.code = code
        self.short_description = short_description
        self.long_description = long_description
        self.no_longer_considered = no_longer_considered
        if not long_description:
            self.long_description = self.short_description

    def __eq__(self, other):
        if isinstance(other, DiscordanceAction):
            return self.code == other.code
        return False

    def __hash__(self):
        return self.code

    def __str__(self):
        return f'({self.code}) {self.short_description}'


DatedDiscordanceAction = typing.NamedTuple('DatedDiscordanceAction', [('action', DiscordanceAction), ('date', Optional[datetime])])


class DiscordanceActionsLog:

    NO_CHANGE = DiscordanceAction(0, "No change in classification")
    NOT_REVIEWED = DiscordanceAction(1, "Not reviewed")

    CHANGE_NO_REASON_GIVEN = DiscordanceAction(10, "Reclassified no reason given")
    CHANGE_UNKNOWN = DiscordanceAction(11, "Reclassified with reason (not supported in current code)")

    CHANGE_DD = DiscordanceAction(20, "Reclassified after discordance discussion")
    CHANGE_ND = DiscordanceAction(21, "Reclassified after summation of data")
    CHANGE_IR = DiscordanceAction(22, "Reclassified after internal review")

    CLINICAL_GROUP_CHANGED = DiscordanceAction(100, "Clinical grouping changed", no_longer_considered=True)
    CLASSIFICATION_WITHDRAWN = DiscordanceAction(101, "Classification withdrawn", no_longer_considered=True)

    ALL_ACTIONS = [
        NO_CHANGE,
        NOT_REVIEWED,
        CHANGE_NO_REASON_GIVEN,
        CHANGE_UNKNOWN,
        CHANGE_DD,
        CHANGE_ND,
        CHANGE_IR,
        CLINICAL_GROUP_CHANGED,
        CLASSIFICATION_WITHDRAWN
    ]

    def __init__(self):
        self.actions = list()
        self.internal_reviewed = None

    def add(self, action: DiscordanceAction, date: Optional[datetime] = None):
        self.actions.append(DatedDiscordanceAction(action=action, date=date))
        self.internal_reviewed = None

    @property
    def no_longer_considered(self):
        for dated_action in self.actions:
            if dated_action.action.no_longer_considered:
                return True
        return False

    def __str__(self):
        return f'internal review = {self.internal_reviewed}, actions = {self.actions}'


class DiscordanceReportClassification(TimeStampedModel):
    report = models.ForeignKey(DiscordanceReport, on_delete=CASCADE)
    classification_original = models.ForeignKey(ClassificationModification, related_name='+', on_delete=CASCADE)
    classification_final = models.ForeignKey(ClassificationModification, related_name='+', on_delete=CASCADE, null=True, blank=True)
    clinical_context_final = models.ForeignKey(ClinicalContext, on_delete=PROTECT, null=True, blank=True)
    withdrawn_final = models.BooleanField(default=False)

    @lazy
    def classification_effective(self) -> ClassificationModification:
        """
        The final state of the classification (or the current if the report is still open)
        """
        if self.classification_final:
            return self.classification_final
        return ClassificationModification.objects.get(
            classification=self.classification_original.classification,
            is_last_published=True
        )

    @lazy
    def clinical_context_effective(self) -> ClinicalContext:
        """
        The final state clinical context of the classification (or the current if the report is still open)
        """
        if self.clinical_context_final:
            return self.clinical_context_final
        return self.classification_original.classification.clinical_context

    @lazy
    def withdrawn_effective(self) -> bool:
        """
        The final "withdrawn" of the classification (or the current if the report is still open)
        """
        if self.classification_final:
            return self.withdrawn_final
        return self.classification_original.classification.withdrawn

    def close(self):
        self.clinical_context_final = self.clinical_context_effective
        self.withdrawn_final = self.withdrawn_effective
        self.classification_final = self.classification_effective
        self.save()

    @lazy
    def action_log(self) -> DiscordanceActionsLog:
        actions = DiscordanceActionsLog()

        if self.withdrawn_effective:
            # if withdrawn, nothing else matters
            actions.add(DiscordanceActionsLog.CLASSIFICATION_WITHDRAWN)
            return actions

        # Need to revisit this, how long do we give for flags to be updated after the discordance report has been closed
        # Currently 1 day, should it be 1 month? (but only if the relevant flag was opened during the discordance period
        # e.g. a significance_change flag opened on day 1, and finally filled in on day 20)
        flag_collection = self.classification_original.classification.flag_collection_safe
        start = self.report.created
        end = self.report.report_completed_date

        # find all flag comments from 1 minute before discordance was recorded
        # (as we want to get specific values)
        # but only for flags that are now closed.
        start -= timedelta(minutes=1)
        if end:
            # add configurable time to get stuff that happened after discordance was completed
            # gives time for people to fill in things
            end += timedelta(days=settings.DISCORDANCE_REPORT_LEEWAY)

        relevant_comments_qs = FlagComment.objects.filter(
            flag__collection=flag_collection,
            # resolution__status=FlagStatus.CLOSED,
            # flag__resolution__status=FlagStatus.CLOSED,
            flag__flag_type__in=[
                classification_flag_types.internal_review,
                classification_flag_types.significance_change,
            ],
        ).order_by('created')

        relevant_comments_qs = relevant_comments_qs.filter(created__gte=start)
        if end:
            relevant_comments_qs = relevant_comments_qs.filter(created__lte=end)

        # note if we had a significant change we treat it as if the record was reviewed
        # e.g. there's no
        had_significant_change = ClinicalSignificance.is_significant_change(
            old_classification=self.classification_original.get(SpecialEKeys.CLINICAL_SIGNIFICANCE),
            new_classification=self.classification_effective.get(SpecialEKeys.CLINICAL_SIGNIFICANCE)
        )

        if self.report.clinical_context != self.clinical_context_effective:
            actions.add(DiscordanceActionsLog.CLINICAL_GROUP_CHANGED)

        # TODO
        # flags = relevant_comments_qs.values_list('flag', flat=True).distinct() # go through in order and overwrite, so latest is most important

        has_reclass_reason = False
        internal_reviewed = False

        processed_flag = set()
        for flag_comment in relevant_comments_qs:  # : :type flag_comment: FlagComment
            # want to update resolution to be
            flag = flag_comment.flag  # : :type flag: Flag
            if flag.id in processed_flag:
                continue
            processed_flag.add(flag.id)

            if flag.flag_type == classification_flag_types.internal_review and flag.resolution.status == FlagStatus.CLOSED:
                internal_reviewed = True  # maybe we should get more info

            elif flag.flag_type == classification_flag_types.significance_change:
                resolution = flag.resolution
                has_reclass_reason = True
                if resolution and resolution.status == FlagStatus.CLOSED:
                    if resolution.id == 'sc_discordance_discussion':
                        actions.add(DiscordanceActionsLog.CHANGE_DD, flag_comment.created)
                    elif resolution.id == 'sc_new_data':
                        actions.add(DiscordanceActionsLog.CHANGE_ND, flag_comment.created)
                    elif resolution.id == 'sc_internal_review':
                        actions.add(DiscordanceActionsLog.CHANGE_IR, flag_comment.created)
                    else:
                        actions.add(DiscordanceActionsLog.CHANGE_UNKNOWN, flag_comment.created)
                else:
                    actions.add(DiscordanceActionsLog.CHANGE_NO_REASON_GIVEN, flag_comment.created)

        if internal_reviewed:
            actions.internal_reviewed = True

        if internal_reviewed and not had_significant_change:
            actions.add(DiscordanceActionsLog.NO_CHANGE)

        if not internal_reviewed and not had_significant_change:
            actions.add(DiscordanceActionsLog.NOT_REVIEWED)

        if had_significant_change and not has_reclass_reason:
            actions.add(DiscordanceActionsLog.CHANGE_NO_REASON_GIVEN)

        return actions
