from collections import defaultdict
from dataclasses import dataclass
from enum import Enum, IntEnum
from functools import total_ordering
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime

from dateutil import relativedelta, tz
from django.http import StreamingHttpResponse
from django.shortcuts import render

from classification.enums import ShareLevel
from classification.models import classification_flag_types, Classification, ClassificationModification
from classification.models.clinical_context_models import CS_TO_NUMBER
from flags.models import FlagComment
from library.utils import delimited_row, IteratableStitcher, IterableTransformer


@total_ordering
@dataclass(frozen=True)
class ClassificationSummary:
    at: datetime
    allele_id: int
    classification_id: int
    clinical_significance: str
    withdrawn: Optional[bool] = None

    def __lt__(self, other):
        return self.at < other.at

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
            withdrawn=self.withdrawn if self.withdrawn is not None else older.withdrawn
        )


class AlleleStatus(IntEnum):
    Empty = 0
    Unique = 1
    Agreement = 2
    Confidence = 3
    Discordant = 4
    Withdrawn = 5


class AlleleSummary:

    def __init__(self):
        self.classifications: Dict[int, ClassificationSummary] = dict()
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

    def withdrawn_count(self):
        return len(self.classifications) - len(self.not_withdrawn)

    def recalculate(self) -> AlleleStatus:
        count = len(self.not_withdrawn)
        if count == 0:
            return AlleleStatus.Empty
        elif count == 1:
            return AlleleStatus.Unique
        else:
            all_values = set()
            all_buckets = set()
            for summary in self.not_withdrawn:
                if bucket := CS_TO_NUMBER.get(summary.clinical_significance, None):
                    all_buckets.add(bucket)
                    all_values.add(summary.clinical_significance)
            if len(all_buckets) >= 2:
                return AlleleStatus.Discordant
            if len(all_values) >= 2:
                return AlleleStatus.Confidence
            return AlleleStatus.Agreement


class ClassificationAccumulationGraph:

    @dataclass
    class _AlleleSummarySnapshot:
        at: datetime
        counts: Dict[AlleleStatus, int]

        @property
        def unique(self):
            return self.counts.get(AlleleStatus.Unique, 0)

        @property
        def agreement(self):
            return self.counts.get(AlleleStatus.Agreement, 0)

        @property
        def confidence(self):
            return self.counts.get(AlleleStatus.Confidence, 0)

        @property
        def discordant(self):
            return self.counts.get(AlleleStatus.Discordant, 0)

        @property
        def withdrawn(self):
            return self.counts.get(AlleleStatus.Withdrawn, 0)

    class _RunningAccumulation:

        def __init__(self):
            self.allele_summaries: Dict[int, AlleleSummary] = defaultdict(AlleleSummary)

        def add_modification(self, summary: ClassificationSummary):
            allele_summary = self.allele_summaries[summary.allele_id]
            allele_summary.add_modification(summary=summary)

        def snapshot(self, at: datetime) -> 'ClassificationAccumulationGraph._AlleleSummarySnapshot':
            counts: Dict[AlleleStatus, int] = defaultdict(int)
            for allele_summary in self.allele_summaries.values():
                status = allele_summary.biggest_status
                if status != AlleleStatus.Empty:
                    counts[status] += allele_summary.classification_count()
                counts[AlleleStatus.Withdrawn] += allele_summary.withdrawn_count()
                allele_summary.reset()
            return ClassificationAccumulationGraph._AlleleSummarySnapshot(at=at, counts=counts)

    @staticmethod
    def withdrawn_iterable():

        flag_collection_id_to_allele_classification: Dict[int, Tuple[int, int]] = dict()

        flag_qs = FlagComment.objects.filter(flag__flag_type=classification_flag_types.classification_withdrawn) \
            .order_by("created") \
            .values_list("resolution__status", "created", "flag__collection_id")

        def transformer(value_tuple):
            status, created, flag_collection_id = value_tuple
            if flag_collection_id not in flag_collection_id_to_allele_classification:
                if classification_match := Classification.objects \
                        .values_list("variant__variantallele__allele_id", "id") \
                        .filter(flag_collection_id=flag_collection_id) \
                        .first():
                    flag_collection_id_to_allele_classification[flag_collection_id] = classification_match
                else:
                    flag_collection_id_to_allele_classification[flag_collection_id] = (0, 0)

            allele_id, classification_id = flag_collection_id_to_allele_classification[flag_collection_id]

            withdrawn = status == 'O'

            return ClassificationSummary(allele_id=allele_id, classification_id=classification_id, at=created,
                                         clinical_significance=None, withdrawn=withdrawn)

        return IterableTransformer(flag_qs, transformer)

    @staticmethod
    def classification_iterable():
        cm_qs_summary = ClassificationModification.objects.filter(published=True, classification__variant__isnull=False,
                                                          share_level__in=ShareLevel.DISCORDANT_LEVEL_KEYS).order_by(
            "created").values_list("classification_id", "created",
                                          "published_evidence__clinical_significance__value",
                                          "classification__variant__variantallele__allele_id")

        def classification_transformer(results_tuple):
            c_id, created, clinical_significance, allele_id = results_tuple
            return ClassificationSummary(allele_id=allele_id, classification_id=c_id, at=created,
                                         clinical_significance=clinical_significance)

        return IterableTransformer(cm_qs_summary, classification_transformer)

    def report(self):

        time_delta = relativedelta.relativedelta(days=7)

        running_accum = self._RunningAccumulation()
        sub_totals: List[ClassificationAccumulationGraph._AlleleSummarySnapshot] = list()

        stitcher = IteratableStitcher(
            iterables=[
                self.withdrawn_iterable(),
                self.classification_iterable()
            ]
        )

        start_date = None
        next_date = None
        for summary in stitcher:
            if not start_date:
                start_date = summary.at.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
                next_date = start_date + time_delta

            while summary.at > next_date:
                sub_totals.append(running_accum.snapshot(at=start_date))
                start_date = next_date
                next_date = next_date + time_delta

            running_accum.add_modification(summary)

        sub_totals.append(running_accum.snapshot(at=start_date))
        sub_totals.append(running_accum.snapshot(at=start_date + time_delta))

        return sub_totals


def download_report(request):
    cag = ClassificationAccumulationGraph()
    report_data = cag.report()

    def iter_report():
        yield delimited_row(["date", "unique", "agreement", "confidence", "discordant", "withdrawn"])
        for row in report_data:
            yield delimited_row([
                row.at.strftime('%Y-%m-%d'),
                row.unique,
                row.agreement,
                row.confidence,
                row.discordant,
                row.withdrawn
            ])

    response = StreamingHttpResponse(iter_report(), content_type="text/csv")
    response['Content-Disposition'] = f'attachment; filename="classification_accumulation_report.csv"'
    return response
