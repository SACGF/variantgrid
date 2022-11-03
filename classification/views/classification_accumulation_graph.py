from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from enum import IntEnum, Enum
from functools import total_ordering
from typing import Dict, List, Any, Optional, Tuple, Set, Union, Iterable

import pandas as pd
from django.http import StreamingHttpResponse

from classification.enums import ShareLevel
from classification.models import classification_flag_types, Classification, ClassificationModification, EvidenceKeyMap
from flags.models import FlagComment
from library.utils import delimited_row, IterableStitcher, IterableTransformer
from snpdb.models import Lab


@total_ordering
@dataclass(frozen=True)
class ClassificationSummary:
    at: datetime
    allele_id: int
    classification_id: int
    clinical_significance: Optional[str]
    withdrawn: Optional[bool] = None
    lab_name: Optional[str] = None

    def __lt__(self, other):
        return self.at < other.at

    @property
    def is_valid(self):
        return self.classification_id != 0

    def merge(self, older: Optional['ClassificationSummary']) -> 'ClassificationSummary':
        if not older:
            return self
        if self.allele_id != older.allele_id or self.classification_id != older.classification_id:
            raise ValueError("Cannot merge summaries for two different classifications")

        return ClassificationSummary(
            at=self.at,
            allele_id=self.allele_id,
            classification_id=self.classification_id,
            clinical_significance=self.clinical_significance or older.clinical_significance,
            withdrawn=self.withdrawn if self.withdrawn is not None else older.withdrawn,
            lab_name=self.lab_name or older.lab_name
        )


class AlleleStatus(IntEnum):
    Empty = 0
    Single = 1
    Agreement = 2
    Confidence = 3
    Discordant = 4
    Withdrawn = 5


allele_status_report = [AlleleStatus.Withdrawn, AlleleStatus.Discordant, AlleleStatus.Confidence, AlleleStatus.Agreement, AlleleStatus.Single]


@dataclass
class AlleleStatusData:
    label: str
    color: str
    line_color: str


allele_status_data = {
    AlleleStatus.Single: AlleleStatusData(label="Single Submitter", color="#ccc", line_color="#888"),
    AlleleStatus.Agreement: AlleleStatusData(label="Agreement", color="#88dd88", line_color="#66bb66"),
    AlleleStatus.Confidence: AlleleStatusData(label="Confidence", color="#fff3cd", line_color="#ddd1ab"),
    AlleleStatus.Discordant: AlleleStatusData(label="Discordant", color="#ebb", line_color="#c99"),
    AlleleStatus.Withdrawn: AlleleStatusData(label="Withdrawn", color="#111", line_color="#000"),
}


class AccumulationReportMode(str, Enum):
    Classification = "classification"
    Allele = "allele"


class AlleleSummary:

    def __init__(self):
        # key is classification ID
        self.classifications: Dict[int, ClassificationSummary] = {}
        self.biggest_status = AlleleStatus.Empty
        self.last_status = AlleleStatus.Empty
        self.not_withdrawn: List[ClassificationSummary] = []

    def add_modification(self, summary: ClassificationSummary):
        self.classifications[summary.classification_id] = summary.merge(
            self.classifications.get(summary.classification_id))
        self.not_withdrawn = [cs for cs in self.classifications.values() if not cs.withdrawn]

        self.last_status = self.recalculate()
        self.biggest_status = max(self.biggest_status, self.last_status)

    def reset(self):
        self.biggest_status = self.last_status

    def classification_count(self) -> int:
        return len(self.not_withdrawn)

    @property
    def multi_lab(self) -> int:
        if len(self.not_withdrawn) <= 1:
            return False
        labs = set()
        for summary in self.not_withdrawn:
            labs.add(summary.lab_name)
            if len(labs) > 1:
                return True
        return False

    def withdrawn_count(self):
        return len(self.classifications) - len(self.not_withdrawn)

    def recalculate(self) -> AlleleStatus:
        count = len(self.not_withdrawn)
        if count == 0:
            return AlleleStatus.Empty
        elif not self.multi_lab:
            return AlleleStatus.Single
        else:
            all_values = set()
            all_buckets = set()
            for summary in self.not_withdrawn:
                if bucket := EvidenceKeyMap.clinical_significance_to_bucket().get(summary.clinical_significance, None):
                    all_buckets.add(bucket)
                    all_values.add(summary.clinical_significance)
            if len(all_buckets) >= 2:
                return AlleleStatus.Discordant
            if len(all_values) >= 2:
                return AlleleStatus.Confidence
            return AlleleStatus.Agreement


class ClassificationAccumulationGraph:

    @dataclass
    class _SummarySnapshot:
        at: datetime
        counts: Dict[Any, int]

    class _RunningAccumulation:

        def __init__(self, mode: AccumulationReportMode):
            # id is allele ID
            self.mode = mode
            self.allele_summaries: Dict[int, AlleleSummary] = defaultdict(AlleleSummary)

        def add_modification(self, summary: ClassificationSummary):
            allele_summary = self.allele_summaries[summary.allele_id]
            allele_summary.add_modification(summary=summary)

        def snapshot(self, at: datetime) -> 'ClassificationAccumulationGraph._SummarySnapshot':
            counts: Dict[Union[AlleleStatus, str], int] = defaultdict(int)
            for allele_summary in self.allele_summaries.values():
                status = allele_summary.biggest_status

                if self.mode == AccumulationReportMode.Classification:
                    if status != AlleleStatus.Empty:
                        counts[status] += allele_summary.classification_count()
                    counts[AlleleStatus.Withdrawn] += allele_summary.withdrawn_count()
                    allele_summary.reset()

                    for summary in allele_summary.not_withdrawn:
                        counts[summary.lab_name] = counts[summary.lab_name] + 1

                elif self.mode == AccumulationReportMode.Allele:
                    if status != AlleleStatus.Empty and allele_summary.classification_count():
                        counts[status] += 1

                    allele_summary.reset()

                    involved_labs: Set[str] = set()
                    for summary in allele_summary.not_withdrawn:
                        involved_labs.add(summary.lab_name)

                    for lab in involved_labs:
                        counts[lab] += 1

            return ClassificationAccumulationGraph._SummarySnapshot(at=at, counts=counts)

    def __init__(self,
                 mode: AccumulationReportMode,
                 shared_only: bool = True,
                 labs: Optional[List[Lab]] = None,
                 time_step: Optional[timedelta] = None
                 ):
        self.mode = mode
        self.shared_only = shared_only
        self.labs = labs
        self.time_delta = time_step or timedelta(days=7)

    @property
    def share_levels(self):
        if self.shared_only:
            return ShareLevel.DISCORDANT_LEVEL_KEYS
        else:
            return ShareLevel.ALL_LEVELS

    def withdrawn_iterable(self) -> IterableTransformer[ClassificationSummary]:
        flag_collection_id_to_allele_classification: Dict[int, Tuple[int, int, Optional[str]]] = {}

        flag_qs = FlagComment.objects.filter(flag__flag_type=classification_flag_types.classification_withdrawn) \
            .order_by("created") \
            .values_list("resolution__status", "created", "flag__collection_id")

        def transformer(value_tuple):
            status, created, flag_collection_id = value_tuple
            if flag_collection_id not in flag_collection_id_to_allele_classification:

                find_class_qs = Classification.objects \
                        .values_list("allele_id", "id", "lab__name") \
                        .filter(flag_collection_id=flag_collection_id) \
                        .filter(share_level__in=self.share_levels)
                if labs := self.labs:
                    find_class_qs = find_class_qs.filter(lab__in=labs)
                else:
                    find_class_qs = find_class_qs.filter(lab__external=False)

                if classification_match := find_class_qs.first():
                    flag_collection_id_to_allele_classification[flag_collection_id] = classification_match
                else:
                    flag_collection_id_to_allele_classification[flag_collection_id] = (0, 0, None)

            allele_id, classification_id, lab_name = flag_collection_id_to_allele_classification[flag_collection_id]

            withdrawn = status == 'O'

            return ClassificationSummary(allele_id=allele_id, classification_id=classification_id, at=created,
                                         clinical_significance=None, withdrawn=withdrawn, lab_name=lab_name)

        return IterableTransformer(flag_qs, transformer)

    def classification_iterable(self):
        cm_qs_summary = ClassificationModification.objects.filter(
            published=True,
            classification__allele__isnull=False,
            share_level__in=self.share_levels).order_by(
            "created").values_list("classification_id", "created",
                                            "published_evidence__clinical_significance__value",
                                            "classification__allele_id",
                                            "classification__lab__name")
        if labs := self.labs:
            cm_qs_summary = cm_qs_summary.filter(classification__lab__in=labs)
        else:
            cm_qs_summary = cm_qs_summary.filter(classification__lab__external=False)

        def classification_transformer(results_tuple):
            c_id, created, clinical_significance, allele_id, lab_name = results_tuple
            if clinical_significance and clinical_significance.startswith('VUS'):
                clinical_significance = 'VUS'

            return ClassificationSummary(allele_id=allele_id, classification_id=c_id, at=created,
                                         clinical_significance=clinical_significance, lab_name=lab_name)

        return IterableTransformer(cm_qs_summary, classification_transformer)

    def report(self) -> List['ClassificationAccumulationGraph._SummarySnapshot']:

        time_delta = self.time_delta
        running_accum = self._RunningAccumulation(mode=self.mode)
        sub_totals: List[ClassificationAccumulationGraph._SummarySnapshot] = []

        stitcher = IterableStitcher[ClassificationSummary](
            iterables=[
                self.withdrawn_iterable(),
                self.classification_iterable()
            ]
        )

        start_date = None
        next_date = None

        summary: ClassificationSummary
        for summary in stitcher:

            if not summary.is_valid:
                continue

            if not start_date:
                start_date = summary.at.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
                next_date = start_date + time_delta

            while summary.at > next_date:
                sub_totals.append(running_accum.snapshot(at=start_date))
                start_date = next_date
                next_date = next_date + time_delta

            running_accum.add_modification(summary)

        if start_date:
            sub_totals.append(running_accum.snapshot(at=start_date))
            sub_totals.append(running_accum.snapshot(at=start_date + time_delta))

        ## if numbers go off, uncomment the following code to test
        # expected_end_classification_ids = set(Classification.objects.filter(allele__isnull=False, lab__in=self.labs, withdrawn=False, share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).values_list('pk', flat=True))
        # for allele_summary in running_accum.allele_summaries.values():
        #     for cs in allele_summary.not_withdrawn:
        #         if cs.classification_id not in expected_end_classification_ids:
        #             raise ValueError(f"{cs.classification_id} was counted but not expected")

        return sub_totals


def _iter_report_list(
        mode: AccumulationReportMode = AccumulationReportMode.Classification,
        shared_only: bool = True,
        labs: Optional[List[Lab]] = None,
        time_step: Optional[timedelta] = None
        ):

    cag = ClassificationAccumulationGraph(mode=mode, shared_only=shared_only, labs=labs, time_step=time_step)
    report_data = cag.report()

    valid_labs = set()
    for row in report_data:
        for key in [key for key in row.counts.keys() if isinstance(key, str)]:
            valid_labs.add(key)
    lab_list = sorted(valid_labs)

    header = ["Date"]
    header.extend([allele_status_data[status].label for status in allele_status_report])
    header.extend(lab_list)

    yield header
    for row in report_data:
        row_date = [row.at.strftime('%Y-%m-%d')]
        for status in allele_status_report:
            row_date.append(row.counts.get(status, 0))
        row_date.extend([row.counts.get(lab, 0) for lab in lab_list])
        yield row_date


def download_report(request):
    mode = AccumulationReportMode.Classification
    if request.GET.get('mode') == 'allele':
        mode = AccumulationReportMode.Allele
    shared_only = True
    if request.GET.get('share') == 'all':
        shared_only = False
    frequency = 7
    if frequency_param := request.GET.get('frequency'):
        frequency = int(frequency_param)

    response = StreamingHttpResponse((delimited_row(r) for r in _iter_report_list(mode=mode, shared_only=shared_only, time_step=timedelta(days=frequency))), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="{mode.value.lower()}_accumulation_report.csv"'
    return response


def get_accumulation_graph_data(mode: AccumulationReportMode = AccumulationReportMode.Classification, labs: Optional[Iterable[Lab]] = None) -> Dict[str, Any]:
    data = list(_iter_report_list(mode=mode, labs=labs))
    header = data[0]
    rows = data[1:]
    df = pd.DataFrame(rows, columns=header)
    labs = df.columns[6:]
    statuses = df.columns[1:5]
    dates = df["Date"].tolist()

    by_lab = []
    by_status = []

    for lab in labs:
        by_lab.append({
            "x": dates,
            "y": df[lab].tolist(),
            "name": lab
        })

    for status in allele_status_report[1:]:
        status_data = allele_status_data[status]
        by_status.append({
            "x": dates,
            "y": df[status_data.label].tolist(),
            "name": status_data.label,
            "fillcolor": status_data.color,
            "line": {
                "color": status_data.line_color
            },
            "stackgroup": "one"
        })

    by_lab.sort(key=lambda t: t["y"][-1], reverse=True)

    return {
        "lab": by_lab,
        "status": by_status
    }
