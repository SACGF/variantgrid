from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import cached_property
from typing import List, Set, Iterable, Optional, Dict

from django.utils import timezone
from frozendict import frozendict
from more_itertools import first

from classification.models import DiscordanceReport, DiscordanceReportTriage, DiscordanceReportTriageStatus, \
    ClassificationModification, ClassificationLabSummary, DiscordanceReportNextStep, DiscordanceReportClassification
from genes.hgvs import CHGVS
from library.django_utils import get_url_from_view_path
from library.utils import ExportRow, export_column, ExportDataType, pretty_label
from snpdb.lab_picker import LabPickerData
from snpdb.models import Lab


class DiscordanceReportRowDataTriagesRowData(ExportRow):

    def __init__(self, discordance_report: DiscordanceReport, perspective: LabPickerData):
        self.discordance_report = discordance_report
        self.perspective = perspective

    @cached_property
    def triages(self) -> List['DiscordanceReportTriage']:
        return list(self.discordance_report.discordancereporttriage_set.order_by('lab__name'))

    def formatted_triages_for_status(self, status: 'DiscordanceReportTriageStatus') -> str:
        triages = [t for t in self.triages if t.triage_status == status]

        def format_triage(t: DiscordanceReportTriage):
            data = f"{t.lab}"
            if note := t.note:
                note = note.replace("\n", "; ")
                data += f": {note}"
            return data

        return "\n".join(format_triage(t) for t in triages)

    @export_column(label="Pending Triage")
    def pending(self):
        return self.formatted_triages_for_status(DiscordanceReportTriageStatus.PENDING)

    @export_column(label="Will Amend")
    def will_amend(self):
        return self.formatted_triages_for_status(DiscordanceReportTriageStatus.REVIEWED_WILL_FIX)

    @export_column(label="No Change")
    def no_change(self):
        return self.formatted_triages_for_status(DiscordanceReportTriageStatus.REVIEWED_SATISFACTORY)

    @export_column(label="Wants Discussion")
    def for_discussion(self):
        return self.formatted_triages_for_status(DiscordanceReportTriageStatus.REVIEWED_WILL_DISCUSS)

    @export_column(label="Low Penetrance etc")
    def low_penetrance_etc(self):
        return self.formatted_triages_for_status(DiscordanceReportTriageStatus.COMPLEX)


class DiscordanceReportRowData(ExportRow):

    def __init__(self, discordance_report: DiscordanceReport, perspective: LabPickerData):
        self.discordance_report = discordance_report
        self.perspective = perspective

    @export_column(label="id")
    def _id(self):
        return self.discordance_report.id

    @export_column(label="Discordance Date", data_type=ExportDataType.date)
    def _discordance_date(self):
        return self.discordance_report.report_started_date

    @export_column(label="URL")
    def _url(self):
        return get_url_from_view_path(self.discordance_report.get_absolute_url())

    @export_column(label="Status")
    def _status(self):
        return self.discordance_report.resolution_text

    @export_column(label="Next Step")
    def _next_step(self):
        if next_step := self.next_step:
            return pretty_label(next_step.name.lower())

    @export_column(label="c.HGVS")
    def _chgvs(self):
        return str(self.c_hgvs)

    @property
    def is_medically_significant(self):
        return self.discordance_report.is_medically_significant

    @export_column(label="Is Medically Significant")
    def _is_medically_significant(self):
        return self.is_medically_significant

    @export_column(label="Lab Significances")
    def _lab_significances(self):
        summaries = self.lab_significances
        rows = []
        for summary in summaries:
            row_parts = [f"{summary.lab}: {summary.clinical_significance_from}"]
            if summary.count > 1:
                row_parts.append(f" x {summary.count}")
            if summary.clinical_significance_from != summary.clinical_significance_to:
                row_parts.append(f" -> {summary.clinical_significance_to}")
                if summary.pending:
                    row_parts.append(" (pending)")
            rows.append(" ".join(row_parts))
        return "\n".join(rows)

    @export_column(label="Conditions")
    def _conditions(self):
        rows = set()
        for cm in self.discordance_report.all_classification_modifications:
            if condition := cm.classification.condition_resolution_obj:
                rows.add(condition.as_plain_text)
        return "\n".join(sorted(rows))

    @export_column(label="Triage", sub_data=DiscordanceReportRowDataTriagesRowData)
    def _triage(self):
        return DiscordanceReportRowDataTriagesRowData(discordance_report=self.discordance_report, perspective=self.perspective)

    @cached_property
    def all_actively_involved_labs(self):
        return self.discordance_report.all_actively_involved_labs

    def all_actively_involved_labs_ids(self) -> str:
        return ";".join(str(lab.pk) for lab in self.all_actively_involved_labs)

    @property
    def is_requiring_attention(self) -> bool:
        return (self.perspective.is_admin_mode or bool(self.all_actively_involved_labs.intersection(self.perspective.your_labs))) \
               and self.discordance_report.resolution is None and not self.discordance_report.is_pending_concordance

    @property
    def is_valid_including_withdraws(self) -> bool:
        return self.perspective.is_admin_mode or bool(set(self.discordance_report.involved_labs.keys()).intersection(self.perspective.your_labs))

    @property
    def other_labs(self) -> Set[Lab]:
        involved_labs = self.all_actively_involved_labs
        if (not self.perspective.has_multi_org_selection) and (selected := self.perspective.selected_labs):
            involved_labs -= selected
        return involved_labs

    @property
    def is_internal(self):
        if self.perspective.has_multi_org_selection:
            return False
        if self.all_actively_involved_labs.issubset(self.perspective.selected_labs):
            return True
        return False

    @property
    def _cm_candidate(self) -> ClassificationModification:
        return first(self._cm_candidates)

    @property
    def _cm_candidates(self) -> Iterable[ClassificationModification]:
        yielded_any = False
        for cm_candidate in self.discordance_report.all_classification_modifications:
            if cm_candidate.classification.lab in self.perspective.your_labs:
                yielded_any = True
                yield cm_candidate
        if not yielded_any:
            yield self.discordance_report.all_classification_modifications[0]

    @property
    def id(self):
        return self.discordance_report.id

    @property
    def date_detected(self) -> datetime:
        return self.discordance_report.report_started_date

    @property
    def date_detected_str(self) -> str:
        date_str = f"{self.date_detected:%Y-%m-%d}"
        if (timezone.now() - self.date_detected) <= timedelta(days=1):
            date_str = f"{date_str} (NEW)"
        return date_str

    @property
    def c_hgvs(self) -> CHGVS:
        return self._cm_candidate.c_hgvs_best(genome_build=self.perspective.genome_build)

    @property
    def c_hgvses(self) -> List[CHGVS]:
        return sorted({candidate.c_hgvs_best(genome_build=self.perspective.genome_build) for candidate in self._cm_candidates})

    @property
    def is_pending_concordance(self):
        return self.discordance_report.is_pending_concordance

    @cached_property
    def lab_significances(self) -> List[ClassificationLabSummary]:
        from classification.models.discordance_lab_summaries import DiscordanceLabSummary
        return DiscordanceLabSummary.for_discordance_report(discordance_report=self.discordance_report, perspective=self.perspective)

    @cached_property
    def next_step(self):
        if not self.is_requiring_attention:
            return None

        is_unanimously_complex = True
        awaiting_your_triage = False
        awaiting_your_amend = False
        awaiting_other_lab = False

        for triage in self.discordance_report.discordancereporttriage_set.select_related('lab').all():
            if not triage.triage_status == DiscordanceReportTriageStatus.COMPLEX:
                is_unanimously_complex = False

            if triage.lab in self.perspective.selected_labs:
                if triage.triage_status == DiscordanceReportTriageStatus.REVIEWED_WILL_FIX:
                    awaiting_your_amend = True
                elif triage.triage_status == DiscordanceReportTriageStatus.PENDING:
                    awaiting_your_triage = True
            else:
                if triage.triage_status in (DiscordanceReportTriageStatus.PENDING, DiscordanceReportTriageStatus.REVIEWED_WILL_FIX):
                    awaiting_other_lab = True

        if is_unanimously_complex:
            return DiscordanceReportNextStep.UNANIMOUSLy_COMPLEX
        elif awaiting_your_triage:
            return DiscordanceReportNextStep.AWAITING_YOUR_TRIAGE
        elif awaiting_your_amend:
            return DiscordanceReportNextStep.AWAITING_YOUR_AMEND
        elif awaiting_other_lab:
            return DiscordanceReportNextStep.AWAITING_OTHER_LAB
        else:
            return DiscordanceReportNextStep.TO_DISCUSS

    def __lt__(self, other):
        def sort_value(x: DiscordanceReportRowData):
            return x.is_medically_significant and x.is_requiring_attention, x.discordance_report.pk
        return sort_value(self) < sort_value(other)


@dataclass(frozen=True)
class DiscordanceReportSummaryCount:
    lab: Optional[Lab]
    count: int

    @property
    def is_internal(self):
        return self.lab is None

    def __lt__(self, other):
        # if we want internal to go first no matter the count
        if self.is_internal != other.is_internal:
            return self.is_internal
        if self.count != other.count:
            return self.count > other.count

        return self.lab < other.lab


class DiscordanceReportTableData:

    def __init__(self, perspective: LabPickerData, summaries: Iterable[DiscordanceReportRowData]):
        self.perspective = perspective
        self.summaries = summaries or []

    def medically_significant_only(self) -> 'DiscordanceReportTableData':
        return DiscordanceReportTableData(
            perspective=self.perspective,
            summaries=[s for s in self.summaries if s.is_medically_significant]
        )

    def __len__(self):
        return len(self.summaries)

    def __bool__(self):
        return bool(self.summaries)

    @property
    def genome_build(self):
        return self.perspective.genome_build

    @cached_property
    def counts(self) -> List[DiscordanceReportSummaryCount]:
        internal_count = 0
        by_lab: Dict[Lab, int] = defaultdict(int)
        for summary in self.summaries:
            if summary.is_internal:
                internal_count += 1
            else:
                for lab in summary.other_labs:
                    by_lab[lab] += 1

        counts = sorted([DiscordanceReportSummaryCount(lab=key, count=value) for key, value in by_lab.items()])
        if internal_count:
            counts.insert(0, DiscordanceReportSummaryCount(lab=None, count=internal_count))
        return counts


@dataclass
class DiscordanceReportCategoriesCounts:
    active: int = 0
    medical: int = 0
    waiting_for_triage: int = 0
    waiting_for_triage_medical: int = 0
    waiting_for_amend: int = 0
    ready_to_discuss: int = 0


class DiscordanceReportCategories:

    def __init__(self, perspective: LabPickerData):
        self.perspective = perspective

        discordant_c = DiscordanceReportClassification.objects \
            .filter(classification_original__classification__lab__in=perspective.selected_labs) \
            .values_list('report_id', flat=True)
        # .filter(classification_original__classification__withdrawn=False)  used to
        self.dr_qs = DiscordanceReport.objects.filter(pk__in=discordant_c)\
            .prefetch_related('discordancereporttriage_set')\
            .prefetch_related('discordancereportclassification_set')\
            .select_related('clinical_context')

    def labs(self) -> List[Lab]:
        return sorted(self.perspective.your_labs)

    def labs_quick_str(self) -> str:
        if len(self.perspective.selected_labs) == 1:
            return str(first(self.perspective.selected_labs))
        else:
            return "your assigned labs"

    @cached_property
    def all_counts(self) -> DiscordanceReportCategoriesCounts:
        counts = DiscordanceReportCategoriesCounts()
        counts.active = len(self.active)
        counts.medical = len([rd for rd in self.active if rd.is_medically_significant])
        triage = self.to_triage_table
        counts.waiting_for_triage = len(triage)
        counts.waiting_for_triage_medical = len(triage.medically_significant_only())
        counts.waiting_for_amend = len(self.to_amend_table)
        counts.ready_to_discuss = len(self.to_discuss_table)
        return counts

    @cached_property
    def unresolved(self) -> List[DiscordanceReportRowData]:
        all_summaries: List[DiscordanceReportRowData] = []
        for dr in self.dr_qs.filter(resolution__isnull=True):
            summary = DiscordanceReportRowData(discordance_report=dr, perspective=self.perspective)
            if summary.is_valid_including_withdraws:
                all_summaries.append(summary)
        all_summaries.sort(reverse=True)
        return all_summaries

    @cached_property
    def active(self) -> List[DiscordanceReportRowData]:
        return [drr for drr in self.unresolved if drr.is_requiring_attention]

    @cached_property
    def historic(self) -> List[DiscordanceReportRowData]:
        historic: List[DiscordanceReportRowData] = []
        # add unresolved but not requiring attention discordance reports (aka pending concordance)
        for summary in self.unresolved:
            if not summary.is_requiring_attention:
                historic.append(summary)

        # add resolved (e.g. concordant/continued discordances)
        for dr in self.dr_qs.filter(resolution__isnull=False):
            summary = DiscordanceReportRowData(discordance_report=dr, perspective=self.perspective)
            if summary.is_valid_including_withdraws:
                historic.append(summary)
        historic.sort(reverse=True)
        return historic

    @cached_property
    def active_by_next_step(self) -> Dict[DiscordanceReportNextStep, List[DiscordanceReportRowData]]:
        by_step = defaultdict(list)
        for row_data in self.active:
            by_step[row_data.next_step].append(row_data)

        return frozendict(by_step)

    @cached_property
    def to_historic_table(self) -> DiscordanceReportTableData:
        return DiscordanceReportTableData(perspective=self.perspective, summaries=self.historic)

    @cached_property
    def to_triage_table(self) -> DiscordanceReportTableData:
        return DiscordanceReportTableData(perspective=self.perspective, summaries=self.active_by_next_step.get(DiscordanceReportNextStep.AWAITING_YOUR_TRIAGE))

    @cached_property
    def to_amend_table(self) -> DiscordanceReportTableData:
        return DiscordanceReportTableData(perspective=self.perspective, summaries=self.active_by_next_step.get(DiscordanceReportNextStep.AWAITING_YOUR_AMEND))

    @cached_property
    def to_waiting_other_table(self) -> DiscordanceReportTableData:
        return DiscordanceReportTableData(perspective=self.perspective, summaries=self.active_by_next_step.get(DiscordanceReportNextStep.AWAITING_OTHER_LAB))

    @cached_property
    def to_complex_table(self) -> DiscordanceReportTableData:
        return DiscordanceReportTableData(perspective=self.perspective, summaries=self.active_by_next_step.get(DiscordanceReportNextStep.UNANIMOUSLy_COMPLEX))

    @cached_property
    def to_discuss_table(self) -> DiscordanceReportTableData:
        return DiscordanceReportTableData(perspective=self.perspective, summaries=self.active_by_next_step.get(DiscordanceReportNextStep.TO_DISCUSS))
