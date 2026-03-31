from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property
from typing import Optional, Iterable, Iterator

from auditlog.context import set_actor
from auditlog.models import LogEntry
from django.contrib.auth.models import User
from django.db.models import QuerySet

from classification.enums import TestingContextBucket, OverlapStatus
from classification.models import ClassificationGrouping, ClassificationResultValue, OverlapContributionStatus, \
    OverlapContribution, OverlapEntrySourceTextChoices, Overlap, OverlapType, OverlapContributionSkew, TriageStatus, \
    TriageNextStep, OverlapContributionLog, OverlapContributionChangeSource
from classification.services.overlap_calculator import calculator_for_value_type, OverlapCalculatorOncPath, \
    OverlapCalculatorClinSig


class OverlapServices:

    @staticmethod
    def update_classification_grouping_overlap_contribution(classification_grouping: ClassificationGrouping):
        if classification_grouping.testing_context in {TestingContextBucket.OTHER, TestingContextBucket.UNKNOWN}:
            # no conflicts for other
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

            # FIXME
            effective_date = classification_grouping.latest_classification_modification.curated_date

            overlap_contribution, created = OverlapContribution.objects.update_or_create(
                source=OverlapEntrySourceTextChoices.CLASSIFICATION,
                allele=classification_grouping.allele,
                classification_grouping=classification_grouping,
                testing_context_bucket=classification_grouping.testing_context,
                tumor_type_category=classification_grouping.tumor_type,
                value_type=value_type,
                defaults={
                    "value": value,
                    "contribution_status": contribution,
                    "effective_date": effective_date
                    # TODO date type
                }
            )
            if created:
                OverlapContributionLog.from_current_overlap_state(overlap_contribution, OverlapContributionChangeSource.UPLOAD).save()
                OverlapServices.link_overlap_contribution(overlap_contribution)
                # make sure this is added to or creates the relevant overlaps
                overlap_contribution.refresh_from_db()
            else:
                last_log = OverlapContributionLog.objects.filter(overlap_contribution=overlap_contribution).order_by('-created').first()
                if value != last_log.value or contribution != last_log.contribution_status or effective_date != last_log.effective_date:
                    OverlapContributionLog.from_current_overlap_state(
                        overlap_contribution,
                        OverlapContributionChangeSource.UPLOAD
                    ).save()

            # now update status of any created overlaps or existing linked overlaps
            for skew in overlap_contribution.overlapcontributionskew_set.select_related('overlap').all():
                OverlapServices.recalc_overlap(skew.overlap)

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
            testing_context=overlap_contribution.testing_context_bucket,
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
            testing_context=None,
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
            status_buckets[skew.contribution.triage_status].append(skew)
            if skew.contribution.triage_status_obj != TriageStatus.NON_INTERACTIVE_THIRD_PARTY:
                all_interactive_skews.append(skew)

        print(status_buckets)
        pending = status_buckets[TriageStatus.PENDING]
        reviewed_will_change = status_buckets[TriageStatus.REVIEWED_WILL_FIX]
        reviewed_will_discuss = status_buckets[TriageStatus.REVIEWED_WILL_DISCUSS]
        reviewed_confident = status_buckets[TriageStatus.REVIEWED_SATISFACTORY]
        reviewed_complex = status_buckets[TriageStatus.COMPLEX]
        non_interactive = status_buckets[TriageStatus.NON_INTERACTIVE_THIRD_PARTY]

        for entry in all_interactive_skews:
            entry.skew_perspective = TriageNextStep.PENDING_CALCULATION

        had_pending_or_changing = False

        # If your triage is pending, your task is always to triage
        if pending:
            had_pending_or_changing = True
            for pend in pending:
                pend.skew_perspective = TriageNextStep.AWAITING_YOUR_TRIAGE
                # if someone else is pending, your task is always to wait on them
                for not_pending in reviewed_will_discuss + reviewed_confident + reviewed_complex:
                    not_pending.skew_perspective = TriageNextStep.AWAITING_OTHER_LAB

        if reviewed_will_change:
            had_pending_or_changing = True
            for change in reviewed_will_change:
                change.skew_perspective = TriageNextStep.AWAITING_YOUR_AMEND
            for not_pending in reviewed_will_discuss + reviewed_confident + reviewed_complex:
                not_pending.skew_perspective = TriageNextStep.AWAITING_OTHER_LAB

        if not had_pending_or_changing:
            # below no one is pending or has an outstanding change, so it's just a matter of working out clashing statuses
            if reviewed_will_discuss:
                # everyone has a status, and at least 1 person said reviewed will discuss, so we're all discussing it now
                for lets_chat in reviewed_will_discuss + reviewed_confident + reviewed_complex:
                    lets_chat.skew_perspective = TriageNextStep.TO_DISCUSS

            elif reviewed_complex or reviewed_confident:
                # no one said discuss, but everyone has a different opinion so time to discuss
                for lets_chat in reviewed_complex + reviewed_confident:
                    lets_chat.skew_perspective = TriageNextStep.TO_DISCUSS

            elif reviewed_complex:
                # No Reviewed Will Discuss
                for complex in reviewed_complex:
                    complex.skew_perspective = TriageNextStep.UNANIMOUSLY_COMPLEX

        # the above should have updated every skew perspective, check below
        for entry in all_interactive_skews:
            if entry.skew_perspective == TriageNextStep.PENDING_CALCULATION:
                raise ValueError("Failed to assign each skew a status")

        OverlapContributionSkew.objects.bulk_update(
            objs=all_interactive_skews,
            fields=['skew_perspective']
        )

    @staticmethod
    def recalc_overlap(overlap: Overlap):
        calculator = calculator_for_value_type(overlap.value_type)
        overlap_status = calculator.calculate_entries(list(overlap.contributions.all()))
        overlap.overlap_status = overlap_status

        if overlap.overlap_type == OverlapType.SINGLE_CONTEXT:
            overlap.valid = True
        elif overlap.overlap_type == OverlapType.CROSS_CONTEXT:
            # cross contexts need at least 2 different contexts to be considered valid
            valid = len(overlap.testing_contexts_objs) > 1
            overlap.valid = valid
        overlap.save()


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

    def __lt__(self, other):
        return self.timestamp < other.timestamp


@dataclass(frozen=True)
class OverlapGrouping:
    single_context_overlap: Optional[Overlap]
    cross_context_overlap: Optional[Overlap]
    value_type: ClassificationResultValue
    user_contribution: OverlapContribution
    contributions: set[OverlapContribution]

    @property
    def user_next_step_single_context(self) -> TriageNextStep:
        if single_context_overlap := self.single_context_overlap:
            if skew := single_context_overlap.overlapcontributionskew_set.filter(contribution=self.user_contribution, overlap__overlap_type=OverlapType.SINGLE_CONTEXT).first():
                return TriageNextStep(skew.skew_perspective)
        return TriageNextStep.PENDING_CALCULATION

    @property
    def user_next_step_cross_context(self) -> TriageNextStep:
        if cross_context_overlap := self.cross_context_overlap:
            if skew := cross_context_overlap.overlapcontributionskew_set.filter(contribution=self.user_contribution, overlap__overlap_type=OverlapType.CROSS_CONTEXT).first():
                return TriageNextStep(skew.skew_perspective)
        return TriageNextStep.PENDING_CALCULATION

    @cached_property
    def change_log(self) -> list[ChangeRow]:
        all_triages = self.contributions | set([self.user_contribution])
        change_rows: list[ChangeRow] = []

        def none_str_to_none(changes: list[str]) -> list[Optional[str]]:
            if changes:
                return [None if value == "None" else value for value in changes]
            return None

        for triage in all_triages:
            triage_log: QuerySet[LogEntry] = LogEntry.objects.get_for_object(triage)
            for entry in triage_log:
                field_changes = []

                new_value_change = none_str_to_none(entry.changes_dict.get("new_value"))
                triage_status_change = none_str_to_none(entry.changes_dict.get("triage_status"))
                if new_value_change or triage_status_change:

                    old_triage_status = TriageStatus.REVIEWED_WILL_FIX
                    new_triage_status = TriageStatus.REVIEWED_WILL_FIX
                    if triage_status_change:
                        old_triage_status = TriageStatus(triage_status_change[0])
                        new_triage_status = TriageStatus(triage_status_change[1])

                    old_new_value = new_value_change[0] if new_value_change else None
                    new_new_value = new_value_change[1] if new_value_change else None

                    old_triage_str = old_triage_status.label
                    if old_new_value:
                        old_new_value = OverlapContribution.pretty_value_for(old_new_value, triage.value_type)
                        old_triage_str += f" ({old_new_value})"

                    new_triage_str = new_triage_status.label
                    if new_new_value:
                        new_new_value = OverlapContribution.pretty_value_for(new_new_value, triage.value_type)
                        new_triage_str += f" ({new_new_value})"

                    field_change = FieldChange("triage_status", old_triage_str, new_triage_str)
                    field_changes.append(field_change)

                for key, value_list in entry.changes_display_dict.items():
                    if key in {"new value", "triage status"}:
                        continue
                    else:
                        print(key)

                    field_change = FieldChange(key, value_list[0], value_list[1])
                    field_changes.append(field_change)

                user: Optional[User] = entry.actor

                change_row = ChangeRow(
                    overlap_contribution=triage,
                    user=user,
                    changes=list(sorted(field_changes)),
                    comment=entry.additional_data.get("comment") if entry.additional_data else None,
                    timestamp=entry.timestamp
                )
                change_rows.append(change_row)
        return list(sorted(change_rows))
        # all_changes = LogEntry.objects.get_for_objects(all_triages)
        # for triage in all_triages:
        #     print(triage.history.latest().changes_dict)
        #
        # print(f"All changes for {all_triages}")
        # for entry in all_changes:
        #     print(entry)
        # print("End changes")

    @staticmethod
    def overlap_grouping_for(classification_grouping: ClassificationGrouping, value_type: ClassificationResultValue, include_cross_context: bool) -> 'OverlapGrouping':
        overlaps = Overlap.objects.filter(
            overlapcontributionskew__contribution__classification_grouping=classification_grouping,
            valid=True,
            value_type=value_type,
            overlap_status__gt=OverlapStatus.NO_CONTRIBUTIONS)
        single_context_overlap = overlaps.filter(overlap_type=OverlapType.SINGLE_CONTEXT).first()
        if not single_context_overlap:
            return None

        contributions = set()
        contributions |= set(single_context_overlap.contributions.filter(contribution_status=OverlapContributionStatus.CONTRIBUTING))

        user_contribution = None
        for contribution in contributions:
            if contribution.classification_grouping_id == classification_grouping.pk:
                user_contribution = contribution
                break

        cross_context_overlap = None
        if include_cross_context:
            cross_context_overlap = overlaps.filter(overlap_type=OverlapType.CROSS_CONTEXT).first()
            if cross_context_overlap:
                contributions |= set(cross_context_overlap.contributions.filter(contribution_status=OverlapContributionStatus.CONTRIBUTING))

        contributions.remove(user_contribution)

        return OverlapGrouping(
            single_context_overlap=single_context_overlap,
            cross_context_overlap=cross_context_overlap,
            value_type=value_type,
            user_contribution=user_contribution,
            contributions=contributions
        )

    @property
    def user_triage(self) -> 'OverlapContribution':
        return self.user_contribution

    @property
    def contribution_count(self):
        return len(self.comparisons) + 1  # the +1 is for the user's contribution

    @cached_property
    def same_context_overlap_status(self) -> OverlapStatus:
        if single_context_overlap := self.single_context_overlap:
            return single_context_overlap.overlap_status_obj
        else:
            return OverlapStatus.NO_CONTRIBUTIONS

    @cached_property
    def cross_context_overlap_status(self) -> OverlapStatus:
        if cross_context_overlap := self.cross_context_overlap:
            return cross_context_overlap.overlap_status_obj
        else:
            return OverlapStatus.NO_CONTRIBUTIONS

    @cached_property
    def comparisons(self) -> list[OverlapEntryCompare]:
        compares = [OverlapEntryCompare(self.user_contribution, other_cont, self.value_type) for other_cont in self.contributions]
        return list(sorted(compares, reverse=True))
