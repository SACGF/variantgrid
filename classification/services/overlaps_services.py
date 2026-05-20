from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property
from json import JSONDecodeError
from typing import Optional, Any
from auditlog.models import LogEntry
from click.decorators import pass_meta_key
from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.utils.timezone import now

from classification.enums import TestingContextBucket, OverlapStatus
from classification.models import ClassificationGrouping, ClassificationResultValue, OverlapContributionStatus, \
    OverlapContribution, OverlapEntrySourceTextChoices, Overlap, OverlapType, OverlapContributionSkew, \
    TriageNextStep, TriageState, EffectiveDate, TriageComment, EffectiveDateType, OverlapDiscordanceNotification
from classification.models.overlaps_enums import TriageStatus
from classification.services.overlap_calculator import calculator_for_value_type, OverlapCalculatorOncPath, \
    OverlapCalculatorClinSig
import json

from classification.signals import send_prepared_discordance_notifications
from snpdb.models import Lab


class OverlapServices:

    @staticmethod
    def update_classification_grouping_overlap_contribution(classification_grouping: ClassificationGrouping):
        if classification_grouping.testing_context in {TestingContextBucket.OTHER, TestingContextBucket.UNKNOWN}:
            # no overlaps for other
            return

        value_types: list[ClassificationResultValue] = [ClassificationResultValue.ONC_PATH]
        if classification_grouping.testing_context != TestingContextBucket.GERMLINE:
            value_types.append(ClassificationResultValue.CLINICAL_SIGNIFICANCE)

        for value_type in value_types:
            calc = calculator_for_value_type(value_type)
            value = calc.value_from_summary(classification_grouping.latest_cached_summary)
            is_comparable = calc.is_comparable_value(value)
            is_shared = classification_grouping.share_level_obj.is_discordant_level

            contribution: OverlapContributionStatus
            if value is None:
                contribution = OverlapContributionStatus.NO_VALUE
            elif not is_shared:
                contribution = OverlapContributionStatus.NOT_SHARED
            elif not is_comparable:
                contribution = OverlapContributionStatus.NON_COMPARABLE_VALUE
            else:
                contribution = OverlapContributionStatus.CONTRIBUTING

            effective_date = EffectiveDate(None, EffectiveDateType.UNKNOWN)
            if lastest_modification := classification_grouping.latest_classification_modification:
                effective_date = lastest_modification.curated_date_check.to_effective_date

            overlap_contribution, created = OverlapContribution.objects.update_or_create(
                source=OverlapEntrySourceTextChoices.CLASSIFICATION,
                allele=classification_grouping.allele,
                classification_grouping=classification_grouping,
                testing_context_bucket=classification_grouping.testing_context,
                tumor_type_category=classification_grouping.tumor_type_category,
                value_type=value_type,
                defaults={
                    "value": value,
                    "contribution_status": contribution,
                    "effective_date": effective_date
                    # TODO date type
                }
            )
            if created:
                OverlapServices.link_overlap_contribution(overlap_contribution)
                # make sure this is added to or creates the relevant overlaps
                overlap_contribution.refresh_from_db()

            # now update status of any created overlaps or existing linked overlaps
            for overlap in overlap_contribution.overlaps:
                OverlapServices.recalc_overlap(overlap)
                OverlapServices.update_skews(overlap)

    @staticmethod
    def link_overlap_contribution(overlap_contribution: OverlapContribution):
        # get single context overlap
        """
        overlap_type = models.TextField(choices=OverlapType.choices)
        value_type = models.TextField(max_length=1, choices=ClassificationResultValue.choices)
        allele = models.ForeignKey(Allele, on_delete=models.CASCADE, null=True, blank=True)  # might be blank for gene symbol wide
        testing_contexts = ArrayField(models.TextField(max_length=1, choices=TestingContextBucket.choices), null=True, blank=True)
        tumor_type_category = models.TextField(null=True, blank=True)  # condition isn't always relevant
        overlap_status = models.IntegerField(choices=OverlapStatus.choices, default=OverlapStatus.NO_CONTRIBUTIONS.value)
        valid = models.BooleanField(default=False)  # if it's cross context but only has contributions from 1 context, or if it's NO_SUBMITTERS it shouldn't be valid

        # have to cache the values
        cached_values = JSONField(null=True, blank=True)
        contributions = models.ManyToManyField(OverlapContribution)
        """

        single_context_overlap, created = Overlap.objects.get_or_create(
            overlap_type=OverlapType.SINGLE_CONTEXT,
            value_type=overlap_contribution.value_type,
            allele=overlap_contribution.allele,
            testing_context_bucket=overlap_contribution.testing_context_bucket,
            tumor_type_category=overlap_contribution.tumor_type_category,
            defaults={
                "overlap_status": OverlapStatus.NO_CONTRIBUTIONS,
                "valid": False
            }
        )

        OverlapContributionSkew.objects.get_or_create(
            overlap=single_context_overlap,
            contribution=overlap_contribution
        )

        cross_context_overlap, created = Overlap.objects.get_or_create(
            overlap_type=OverlapType.CROSS_CONTEXT,
            value_type=overlap_contribution.value_type,
            allele=overlap_contribution.allele,
            testing_context_bucket=None,
            tumor_type_category=None,
            defaults={
                "overlap_status": OverlapStatus.NO_CONTRIBUTIONS,
                "valid": False
            }
        )

        OverlapContributionSkew.objects.get_or_create(
            overlap=cross_context_overlap,
            contribution=overlap_contribution
        )

    @staticmethod
    def update_skews(overlap: Overlap):
        status_buckets: defaultdict[TriageStatus, list[OverlapContributionSkew]] = defaultdict(list)
        all_interactive_skews: list[OverlapContributionSkew] = []
        for skew in overlap.overlapcontributionskew_set.all():
            # move skews into - user has done something, user is waiting on something
            status_buckets[skew.contribution.triage_state.status].append(skew)
            if skew.contribution.triage_state.status != TriageStatus.NON_INTERACTIVE_THIRD_PARTY:
                all_interactive_skews.append(skew)

        pending = status_buckets[TriageStatus.PENDING]
        reviewed_will_change = status_buckets[TriageStatus.REVIEWED_WILL_FIX]
        reviewed_will_discuss = status_buckets[TriageStatus.REVIEWED_WILL_DISCUSS]
        reviewed_confident = status_buckets[TriageStatus.REVIEWED_SATISFACTORY]
        reviewed_complex = status_buckets[TriageStatus.COMPLEX]
        non_interactive = status_buckets[TriageStatus.NON_INTERACTIVE_THIRD_PARTY]

        for entry in all_interactive_skews:
            entry.next_step = TriageNextStep.PENDING_CALCULATION

        had_pending_or_changing = False

        # If your triage is pending, your task is always to triage
        if pending:
            had_pending_or_changing = True

            partially_triaged = False
            # now while there are pending, see if there are records that aren't pending
            for not_pending in reviewed_will_discuss + reviewed_confident + reviewed_complex:
                not_pending.next_step = TriageNextStep.AWAITING_OTHER_LAB
                partially_triaged = True

            for pend in pending:
                if partially_triaged:
                    pend.next_step = TriageNextStep.AWAITING_YOUR_TRIAGE_OTHERS_TRIAGED
                else:
                    pend.next_step = TriageNextStep.AWAITING_YOUR_TRIAGE

        if reviewed_will_change:
            had_pending_or_changing = True
            for change in reviewed_will_change:
                change.next_step = TriageNextStep.AWAITING_YOUR_AMEND
            for not_pending in reviewed_will_discuss + reviewed_confident + reviewed_complex:
                not_pending.next_step = TriageNextStep.AWAITING_OTHER_LAB

        if not had_pending_or_changing:
            # below no one is pending or has an outstanding change, so it's just a matter of working out clashing statuses
            if reviewed_will_discuss:
                # everyone has a status, and at least 1 person said reviewed will discuss, so we're all discussing it now
                for lets_chat in reviewed_will_discuss + reviewed_confident + reviewed_complex:
                    lets_chat.next_step = TriageNextStep.TO_DISCUSS

            elif reviewed_complex or reviewed_confident:
                # no one said discuss, but everyone has a different opinion so time to discuss
                for lets_chat in reviewed_complex + reviewed_confident:
                    lets_chat.next_step = TriageNextStep.TO_DISCUSS

            elif reviewed_complex:
                # No Reviewed Will Discuss
                for complex in reviewed_complex:
                    complex.next_step = TriageNextStep.UNANIMOUSLY_COMPLEX

        # the above should have updated every skew perspective, check below
        for entry in all_interactive_skews:
            if entry.next_step == TriageNextStep.PENDING_CALCULATION:
                raise ValueError("Failed to assign each skew a status")

        OverlapContributionSkew.objects.bulk_update(
            objs=all_interactive_skews,
            fields=['next_step']
        )

    @staticmethod
    def recalc_overlap(overlap: Overlap):
        calculator = calculator_for_value_type(overlap.value_type)
        overlap_status = calculator.calculate_entries(list(overlap.contributions.all()))
        old_overlap_status = overlap.overlap_status
        overlap.overlap_status = overlap_status
        overlap_changed = False
        if overlap_status != old_overlap_status:
            overlap.overlap_status_change_timestamp = now()
            overlap_changed = True
        if overlap.overlap_type == OverlapType.SINGLE_CONTEXT:
            overlap.valid = True
        elif overlap.overlap_type == OverlapType.CROSS_CONTEXT:
            # cross contexts need at least 2 different contexts to be considered valid
            valid = len(overlap.testing_contexts_objs) > 1
            overlap.valid = valid
        overlap.save()
        if overlap_changed:
            OverlapServices.overlap_status_changed(overlap, old_overlap_status)

    @staticmethod
    def overlap_status_changed(overlap: Overlap, old_status: OverlapStatus):
        if not settings.DISCORDANCE_ENABLED:
            return

        if overlap.overlap_type == OverlapType.CROSS_CONTEXT:
            return  # don't notify when cross context become discordant

        new_status = overlap.overlap_status

        # see if it's worth notifying anyone
        if not (new_status.is_discordant ^ old_status.is_discordant):
            # only care if we're going from discordant to not discordant or vice versa
            # an overlap could change from concordant to discordant back to concordant
            # but we just check if the notification is worth sending
            return

        odn, created = OverlapDiscordanceNotification.objects.get_or_create(
            overlap=overlap,
            notification_sent_date=None,
            defaults={"old_status": old_status, "new_status": new_status}
        )
        if not created:
            odn.new_status = new_status
            odn.save()

        # this will send emails right away if we're no in an import, otherwise at the end of an import
        send_prepared_discordance_notifications()


@dataclass
class OverlapEntryCompare:
    entry_1: OverlapContribution
    entry_2: OverlapContribution
    value_type: ClassificationResultValue

    @cached_property
    def comparison(self) -> OverlapStatus:
        if self.value_type == ClassificationResultValue.ONC_PATH:
            return OverlapCalculatorOncPath.calculate_entries([self.entry_1, self.entry_2])
        elif self.value_type == ClassificationResultValue.CLINICAL_SIGNIFICANCE:
            return OverlapCalculatorClinSig.calculate_entries([self.entry_1, self.entry_2])
        else:
            raise ValueError(f"Unsupported value type {self.value_type}")

    @property
    def is_cross_context(self) -> bool:
        return self.entry_1.testing_context_full != self.entry_2.testing_context_full

    @property
    def sort_value(self):
        return (
            not self.is_cross_context,
            self.comparison,
            self.entry_2.testing_context_bucket
        )

    def __lt__(self, other):
        return self.sort_value < other.sort_value

    @property
    def triage(self) -> OverlapContribution:
        return self.entry_2


@dataclass
class FieldChange:
    field: str
    old_value: str
    new_value: str

    def __lt__(self, other):
        return self.field < other.field


@dataclass
class ChangeRow:
    overlap_contribution: OverlapContribution
    user: Optional[User]
    changes: list[FieldChange]
    comment: Optional[str]
    timestamp: datetime
    is_new_record: bool

    def __lt__(self, other):
        return self.timestamp < other.timestamp


@dataclass(frozen=True)
class OverlapContributionPerspective:
    overlap_contribution: OverlapContribution
    is_cross_context: bool = False
    is_user_lab: bool = False

    @property
    def _sort_index(self):
        return self.is_cross_context, self.is_user_lab, self.overlap_contribution

    def __lt__(self, other):
        return self._sort_index < other._sort_index


@dataclass(frozen=True)
class OverlapGrouping3:
    overlap: Overlap
    user: User

    @cached_property
    def rows(self) -> list[OverlapContributionPerspective]:
        perspectives = []

        user_labs = Lab.valid_labs_qs(self.user)
        for overlap_contribution in self.overlap.contributions:
            lab = None
            if cg := overlap_contribution.classification_grouping:
                lab = cg.lab

            perspectives.append(OverlapContributionPerspective(
                overlap_contribution=overlap_contribution,
                is_cross_context=False,
                is_user_lab=not self.user.is_superuser and lab in user_labs
            ))

        if allele_id := self.overlap.allele_id:
            if cross_context := Overlap.objects.filter(allele_id=allele_id, overlap_type=OverlapType.CROSS_CONTEXT, valid=True).first():
                cross_context_contributions = set(cross_context.contributions) - set(self.overlap.contributions_list)
                for overlap_contribution in list(sorted(cross_context_contributions)):
                    perspectives.append(OverlapContributionPerspective(
                        overlap_contribution=overlap_contribution,
                        is_cross_context=True,
                        is_user_lab=not self.user.is_superuser and overlap_contribution.lab in user_labs
                    ))

        return list(sorted(perspectives))

    @property
    def value_type(self) -> ClassificationResultValue:
        return self.overlap.value_type

    @staticmethod
    def tidy_change(overlap_contribution: OverlapContribution, field_name: str, value: Any):
        # just used for the changelog

        def to_json(thing: Any):
            if isinstance(thing, str):
                # sometimes json is double encoded in audit
                thing = json.loads(thing)
                if isinstance(thing, str):
                    thing = json.loads(thing)
            return thing

        if value == "None":
            return None
        try:
            match field_name:
                case "triage_state":
                    return TriageState.from_dict(to_json(value))
                case "effective_date":
                    return EffectiveDate.from_dict(to_json(value))
                case "comment":
                    return TriageComment.from_dict(to_json(value))
        except JSONDecodeError as jer:
            return f"Decoding Error: {value}"
        return value

    @cached_property
    def change_log(self) -> list[ChangeRow]:
        all_triages = self.overlap.contributions_list
        change_rows: list[ChangeRow] = []

        for triage in all_triages:
            triage_log: QuerySet[LogEntry] = LogEntry.objects.get_for_object(triage).order_by('timestamp')

            buffer: list[LogEntry] = []
            for entry in triage_log:
                is_new_record = False
                if (id_change := entry.changes_dict.get("id")) and id_change[0] == 'None':
                    is_new_record = True

                    contributing_status = entry.changes_dict.get("contribution_status")
                    if not contributing_status or contributing_status[1] not in (OverlapContributionStatus.CONTRIBUTING, OverlapContributionStatus.NON_COMPARABLE_VALUE):
                        buffer.append(entry)
                        continue
                elif buffer:
                    merge_buffer = False
                    if contributing_status := entry.changes_dict.get("contribution_status"):
                        if contributing_status[1] in (OverlapContributionStatus.CONTRIBUTING, OverlapContributionStatus.NON_COMPARABLE_VALUE):
                            merge_buffer = True

                    if not merge_buffer:
                        buffer.append(entry)

                latest_values = {}
                if buffer:
                    is_new_record = True  # as in this is the first time we're displaying the record
                    for buffered_entry in reversed([entry] + buffer):
                        for key, value_list in buffered_entry.changes_dict.items():
                            if key not in latest_values:
                                latest_values[key] = value_list
                    buffer.clear()
                else:
                    latest_values = entry.changes_dict

                comment: Optional[TriageComment] = None
                field_changes = []
                for key, value_list in latest_values.items():

                    if key not in ("effective_date", "triage_status", "value", "triage_state"):
                        continue
                    if is_new_record and key == "triage_state":
                        continue  # triage_state should always start as pending

                    old_value = OverlapGrouping3.tidy_change(triage, key, value_list[0])
                    new_value = OverlapGrouping3.tidy_change(triage, key, value_list[1])

                    if key == "comment":
                        comment = new_value
                    else:
                        field_change = FieldChange(key, old_value, new_value)
                        field_changes.append(field_change)

                if field_changes:
                    user: Optional[User] = entry.actor
                    change_row = ChangeRow(
                        overlap_contribution=triage,
                        user=user,
                        changes=list(sorted(field_changes)),
                        comment=comment,
                        timestamp=entry.timestamp,
                        is_new_record=is_new_record
                    )
                    change_rows.append(change_row)

        return list(sorted(change_rows))
