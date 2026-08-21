import statistics
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from math import ceil, log10
from functools import cached_property
from typing import Any, Optional

from dateutil.relativedelta import relativedelta

from django.core.cache import cache
from django.db import connection
from django.db.models import Count, Min, OuterRef, QuerySet, Subquery
from django.db.models.functions import ExtractYear
from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render
from django.utils.timezone import localdate

from classification.enums import (ClinicalSignificance, CriteriaEvaluation, LabExternalFilter,
                                 ReclassificationEventType, SpecialEKeys)
from classification.models import (
    Classification,
    ClassificationModification,
    EvidenceKey,
    EvidenceKeyMap,
    ReclassificationBuildResult,
    ReclassificationEvent,
    ReclassificationEventBuilder,
    classification_flag_types,
)
from classification.tasks.classification_reclassification_tasks import reclassification_events_update
from flags.models import Flag
from library.django_utils import require_superuser
from snpdb.models import Lab, Organization

# top to bottom / left to right on every chart and table
SIGNIFICANCE_ORDER = [
    ClinicalSignificance.PATHOGENIC,
    ClinicalSignificance.LIKELY_PATHOGENIC,
    ClinicalSignificance.VUS,
    ClinicalSignificance.LIKELY_BENIGN,
    ClinicalSignificance.BENIGN,
]

BENIGN_DIRECTION_COLOUR = ClinicalSignificance.chart_colour(ClinicalSignificance.BENIGN)
PATHOGENIC_DIRECTION_COLOUR = ClinicalSignificance.chart_colour(ClinicalSignificance.PATHOGENIC)
# counts that carry no benign/pathogenic meaning stay outside the significance language
NEUTRAL_COUNT_COLOUR = "#6c757d"
# the matrix shades by how busy a cell is, which is a different thing from significance, so it gets
# its own single hue rather than borrowing the pill colours
MATRIX_HEAT_COLOUR = "#4a6f8a"

RECLASSIFICATION_PAGE_BUILD_LIMIT = 5000

TOP_GENE_COUNT = 30
EVIDENCE_MOVEMENT_CACHE_SECONDS = 60 * 60
# the significance itself changes by definition, the dates get stamped alongside it, and the report and
# internal notes keys track a record's paperwork rather than the evidence behind the call
EVIDENCE_KEYS_EXCLUDED = {
    SpecialEKeys.CLINICAL_SIGNIFICANCE,
    SpecialEKeys.CURATION_DATE,
    SpecialEKeys.CURATION_VERIFIED_DATE,
    "internal_use",
    "report_date",
    "report_id",
    "report_type",
}
# criteria always stand on their own, the busiest few other keys join them, and the tail folds into one row
TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT = 8
OTHER_EVIDENCE_GROUP = "__other__"
EVIDENCE_APPLIED_COLOUR = "#4a8f5b"
EVIDENCE_UNAPPLIED_COLOUR = "#b5564e"
EVIDENCE_STRENGTHENED_COLOUR = "#84b590"
EVIDENCE_WEAKENED_COLOUR = "#d59a95"
EVIDENCE_CHANGED_COLOUR = "#4f7ca6"

# curation effort and calls moving are counted side by side all over the page, and neither of them is
# a benign or pathogenic direction, so they stay out of the significance palette
REVIEW_COLOUR = "#5b8f8a"
CHANGE_COLOUR = "#a8743f"

SURVIVAL_ORIGIN_YEARS_BACK = 5
DAYS_PER_YEAR = 365.25
# most records that move do it inside the first couple of years, so the life table steps every 6 months
SURVIVAL_INTERVAL_YEARS = 0.5
SERIAL_YEAR_CAP = 4

# the days axis is logarithmic - most records move within a few months, and a year either side of the
# two year mark is not a distinction anybody reads off a chart
LOG_BINS_PER_DECADE = 8
DAY_TICKS = [(0, "same day"), (7, "1w"), (30, "1m"), (90, "3m"), (182, "6m"), (365, "1y"),
             (730, "2y"), (1826, "5y"), (3652, "10y")]


def significance_label(clinical_significance: Optional[str]) -> str:
    return ClinicalSignificance.SHORT_LABELS.get(clinical_significance) or "?"


@dataclass(frozen=True)
class Distribution:
    """ The shape of a set of day counts, for the time to reclassification summary table """
    label: str
    clinical_significance: str
    days: list[int]

    @property
    def count(self) -> int:
        return len(self.days)

    @property
    def median(self) -> Optional[int]:
        return round(statistics.median(self.days)) if self.days else None

    @property
    def p25(self) -> Optional[int]:
        return self._quantile(0.25)

    @property
    def p75(self) -> Optional[int]:
        return self._quantile(0.75)

    @property
    def longest(self) -> Optional[int]:
        return max(self.days) if self.days else None

    def _quantile(self, fraction: float) -> Optional[int]:
        if not self.days:
            return None
        ordered = sorted(self.days)
        return ordered[min(int(fraction * len(ordered)), len(ordered) - 1)]

    @property
    def css_class(self) -> str:
        return ClinicalSignificance.css_class(self.clinical_significance)

    @property
    def colour(self) -> str:
        return ClinicalSignificance.chart_colour(self.clinical_significance)


@dataclass(frozen=True)
class SignificanceTransition:
    from_clinical_significance: str
    to_clinical_significance: str
    count: int


@dataclass(frozen=True)
class SignificanceMatrixCell:
    count: int
    is_self: bool
    background: str


@dataclass(frozen=True)
class SignificanceMatrixRow:
    clinical_significance: str
    cells: list[SignificanceMatrixCell]

    @property
    def label(self) -> str:
        return significance_label(self.clinical_significance)

    @property
    def css_class(self) -> str:
        return ClinicalSignificance.css_class(self.clinical_significance)

    @property
    def total(self) -> int:
        return sum(cell.count for cell in self.cells)


@dataclass(frozen=True)
class YearlyActivity:
    """ One year of curation effort against the catalogue that was there to be revisited """
    year: int
    population: int
    reviewed: int
    reclassified: int

    @property
    def reviewed_percent(self) -> Optional[float]:
        return round(100 * self.reviewed / self.population, 2) if self.population else None

    @property
    def reclassified_percent(self) -> Optional[float]:
        return round(100 * self.reclassified / self.population, 2) if self.population else None

    @property
    def reviews_per_change(self) -> Optional[float]:
        """ How many records a curator looked at for each one whose call moved """
        return round(self.reviewed / self.reclassified, 1) if self.reclassified else None


@dataclass(frozen=True)
class SerialCount:
    """ How many records were touched in one year, two years, three, and so on """
    label: str
    reviewed: int
    reclassified: int


@dataclass(frozen=True)
class SurvivalCurve:
    """
    Life table for one kind of event over a fixed cohort, on a yearly grid. Records leaving observation -
    withdrawn, or the window running out - are censored rather than counted as never having the event.
    """
    label: str
    event_label: str
    colour: str
    survival: list[float]
    """ Share of the cohort still waiting for the event, at each year of the grid, starting at 1.0 """
    events: int
    censored: int

    @property
    def half_life_years(self) -> Optional[float]:
        for interval, remaining in enumerate(self.survival):
            if remaining <= 0.5:
                return round(interval * SURVIVAL_INTERVAL_YEARS, 1)
        return None

    @property
    def ever_percent(self) -> float:
        return round(100 * (1 - self.survival[-1]), 1) if self.survival else 0.0


@dataclass(frozen=True)
class PointsTransition:
    """ The ACMG points a set of records travelled to make the same move between two significances """
    from_clinical_significance: str
    to_clinical_significance: str
    deltas: list[int]

    @property
    def label(self) -> str:
        return f"{significance_label(self.from_clinical_significance)} \u2192 " \
               f"{significance_label(self.to_clinical_significance)}"

    @property
    def count(self) -> int:
        return len(self.deltas)

    @property
    def median(self) -> Optional[float]:
        return round(statistics.median(self.deltas), 1) if self.deltas else None

    @property
    def colour(self) -> str:
        moved_benign = ClinicalSignificance.distance(
            self.from_clinical_significance, self.to_clinical_significance) or 0
        return BENIGN_DIRECTION_COLOUR if moved_benign > 0 else PATHOGENIC_DIRECTION_COLOUR


@dataclass(frozen=True)
class SurvivalAnalysis:
    """
    A cohort taken on one date, followed to the end of the window. The gap between the two curves is the
    reanalysis that left the call where it was.
    """
    cohort_size: int
    origin: date
    window_end: date
    origin_significance: str
    intervals: list[float]
    """ Years since the origin at each life table step """
    reviewed: SurvivalCurve
    reclassified: SurvivalCurve

    @property
    def significance_label(self) -> str:
        return significance_label(self.origin_significance)

    @property
    def css_class(self) -> str:
        return ClinicalSignificance.css_class(self.origin_significance)

    @property
    def span_years(self) -> float:
        return self.intervals[-1] if self.intervals else 0


@dataclass(frozen=True)
class LabActivity:
    """ One lab's maintenance of its own back catalogue, for the league table """
    lab_name: str
    organization_name: str
    held: int
    reviewed: int
    reclassified: int
    baseline: float
    """ Reclassifications a lab of this size would make at the median per record rate """
    is_selected: bool

    @property
    def reviewed_percent(self) -> Optional[float]:
        return round(100 * self.reviewed / self.held, 2) if self.held else None

    @property
    def reclassified_percent(self) -> Optional[float]:
        return round(100 * self.reclassified / self.held, 2) if self.held else None

    @property
    def reviews_per_change(self) -> Optional[float]:
        return round(self.reviewed / self.reclassified, 1) if self.reclassified else None

    @property
    def residual(self) -> float:
        return round(self.reclassified - self.baseline, 1)


@dataclass(frozen=True)
class GeneBurden:
    gene_symbol: str
    vus_count: int
    total_count: int

    @property
    def percent(self) -> float:
        return round(100 * self.vus_count / self.total_count, 1) if self.total_count else 0.0


@dataclass(frozen=True)
class EvidenceMovement:
    """
    How often one evidence key moved across a set of reclassifications. A criterion coming on, coming off,
    or shifting strength is counted as evidence moving, everything else is a plain change.
    """
    key: str
    label: str
    applied: int
    unapplied: int
    strengthened: int
    weakened: int
    changed: int
    is_criteria: bool

    @property
    def total(self) -> int:
        return self.applied + self.unapplied + self.strengthened + self.weakened + self.changed


class ReclassificationAnalytics:
    """
    Every chart on /classification/reclassification_analytics reads from here, so they all share the
    organisation / lab / date filters the user picked. @see ReclassificationEvent
    """

    def __init__(self, organization: Optional[Organization], lab: Optional[Lab],
                 date_from: Optional[date], date_to: Optional[date], lab_external: LabExternalFilter,
                 origin: Optional[date], origin_significance: str):
        self.organization = organization
        self.lab = lab
        self.date_from = date_from
        self.date_to = date_to
        self.lab_external = lab_external
        self._origin = origin
        self.origin_significance = origin_significance

    @staticmethod
    def from_request(request: HttpRequest) -> 'ReclassificationAnalytics':
        def get_date(param: str) -> Optional[date]:
            if raw := request.GET.get(param):
                try:
                    return date.fromisoformat(raw)
                except ValueError:
                    pass
            return None

        organization = None
        lab = None
        if organization_id := request.GET.get('organization'):
            organization = Organization.objects.filter(pk=organization_id).first()
        if lab_id := request.GET.get('lab'):
            lab = Lab.objects.filter(pk=lab_id).first()
            if lab:
                organization = lab.organization

        lab_external = LabExternalFilter.ALL
        if raw_lab_external := request.GET.get('lab_external'):
            if raw_lab_external in LabExternalFilter.values:
                lab_external = LabExternalFilter(raw_lab_external)

        origin_significance = ClinicalSignificance.VUS
        if raw_significance := request.GET.get('origin_significance'):
            if raw_significance in ClinicalSignificance.LABELS:
                origin_significance = raw_significance

        return ReclassificationAnalytics(
            organization=organization,
            lab=lab,
            date_from=get_date('date_from'),
            date_to=get_date('date_to'),
            lab_external=lab_external,
            origin=get_date('origin'),
            origin_significance=origin_significance)

    # -- querysets every chart builds on --------------------------------------------------------

    @cached_property
    def timeline_qs(self) -> QuerySet[ReclassificationEvent]:
        """ The whole timeline for the selected labs, unbounded by date - denominators need the history """
        qs = ReclassificationEvent.objects.all()
        if self.lab:
            qs = qs.filter(lab=self.lab)
        elif self.organization:
            qs = qs.filter(lab__organization=self.organization)
        if lab_external_q := self.lab_external.filter_q("lab"):
            qs = qs.filter(lab_external_q)
        return qs

    @cached_property
    def events_qs(self) -> QuerySet[ReclassificationEvent]:
        qs = self.timeline_qs
        if self.date_from:
            qs = qs.filter(reclassified_date__gte=self.date_from)
        if self.date_to:
            qs = qs.filter(reclassified_date__lte=self.date_to)
        return qs

    @cached_property
    def reclassifications_qs(self) -> QuerySet[ReclassificationEvent]:
        return self.events_qs.filter(event_type=ReclassificationEventType.RECLASSIFICATION)

    @cached_property
    def reviews_qs(self) -> QuerySet[ReclassificationEvent]:
        """ Every time a curator revisited a record, whether or not the call moved """
        return self.events_qs.filter(event_type__in=ReclassificationEventType.review_types())

    @cached_property
    def reclassification_count(self) -> int:
        return self.reclassifications_qs.count()

    @cached_property
    def review_count(self) -> int:
        return self.reviews_qs.count()

    @cached_property
    def reevaluation_count(self) -> int:
        return self.events_qs.filter(event_type=ReclassificationEventType.REEVALUATION).count()

    @cached_property
    def window_end(self) -> date:
        return self.date_to or localdate()

    @cached_property
    def origin(self) -> date:
        """ Where the survival cohort is taken from - the start of the date filter, or five years back """
        return self._origin or self.date_from \
            or (self.window_end - relativedelta(years=SURVIVAL_ORIGIN_YEARS_BACK))

    @cached_property
    def classification_count(self) -> int:
        return self.timeline_qs.values('classification_id').distinct().count()

    @cached_property
    def labs(self) -> QuerySet[Lab]:
        qs = Lab.objects.filter(organization=self.organization) if self.organization else Lab.objects.all()
        return qs.order_by('organization__name', 'name')

    @cached_property
    def selected_labs(self) -> list[Lab]:
        """ The labs the organisation / lab / internal filters leave holding a timeline """
        lab_ids = self.timeline_qs.values_list('lab', flat=True).distinct()
        return list(Lab.objects.filter(pk__in=lab_ids).select_related('organization')
                    .order_by('organization__name', 'name'))

    @cached_property
    def organizations(self) -> QuerySet[Organization]:
        return Organization.objects.filter(active=True).order_by('name')

    # -- 1. Sankey ------------------------------------------------------------------------------

    @cached_property
    def significance_transitions(self) -> list[SignificanceTransition]:
        counts = self.reclassifications_qs \
            .values('from_clinical_significance', 'to_clinical_significance') \
            .annotate(total=Count('pk')).order_by('-total')
        return [SignificanceTransition(from_clinical_significance=row['from_clinical_significance'],
                                       to_clinical_significance=row['to_clinical_significance'],
                                       count=row['total']) for row in counts]

    @cached_property
    def significance_flow(self) -> dict[str, Any]:
        # a source column then a target column, so a significance that both loses and gains records
        # shows up on each side rather than looping back on itself
        source_index = {cs: index for index, cs in enumerate(SIGNIFICANCE_ORDER)}
        target_index = {cs: index + len(SIGNIFICANCE_ORDER) for index, cs in enumerate(SIGNIFICANCE_ORDER)}

        sources, targets, values, link_colours = [], [], [], []
        source_totals = defaultdict(int)
        target_totals = defaultdict(int)
        for transition in self.significance_transitions:
            from_cs = transition.from_clinical_significance
            to_cs = transition.to_clinical_significance
            if from_cs not in source_index or to_cs not in target_index:
                continue
            sources.append(source_index[from_cs])
            targets.append(target_index[to_cs])
            values.append(transition.count)
            source_totals[from_cs] += transition.count
            target_totals[to_cs] += transition.count
            moved_benign = ClinicalSignificance.distance(from_cs, to_cs) or 0
            link_colours.append(self._translucent(
                BENIGN_DIRECTION_COLOUR if moved_benign > 0 else PATHOGENIC_DIRECTION_COLOUR))

        return {
            "labels": [f"{significance_label(cs)} from" for cs in SIGNIFICANCE_ORDER] +
                      [f"{significance_label(cs)} to" for cs in SIGNIFICANCE_ORDER],
            "node_colours": [ClinicalSignificance.chart_colour(cs) for cs in SIGNIFICANCE_ORDER] * 2,
            # the template stacks these into node positions, so both columns run pathogenic to benign
            "node_totals": [source_totals[cs] for cs in SIGNIFICANCE_ORDER] +
                           [target_totals[cs] for cs in SIGNIFICANCE_ORDER],
            "sources": sources,
            "targets": targets,
            "values": values,
            "link_colours": link_colours,
        }

    @cached_property
    def significance_matrix(self) -> list[SignificanceMatrixRow]:
        counts = {(transition.from_clinical_significance, transition.to_clinical_significance): transition.count
                  for transition in self.significance_transitions}
        busiest = max([counts.get((from_cs, to_cs), 0)
                       for from_cs in SIGNIFICANCE_ORDER for to_cs in SIGNIFICANCE_ORDER] or [0])
        return [
            SignificanceMatrixRow(
                clinical_significance=from_cs,
                cells=[SignificanceMatrixCell(
                    count=counts.get((from_cs, to_cs), 0),
                    is_self=from_cs == to_cs,
                    background=self._heat(counts.get((from_cs, to_cs), 0), busiest))
                    for to_cs in SIGNIFICANCE_ORDER])
            for from_cs in SIGNIFICANCE_ORDER
        ]

    @staticmethod
    def _heat(count: int, busiest: int) -> str:
        """ Square rooted so the quiet cells still read, and capped well short of the pill colours """
        if not count or not busiest:
            return ""
        alpha = round(0.06 + 0.44 * (count / busiest) ** 0.5, 3)
        return ReclassificationAnalytics._translucent(MATRIX_HEAT_COLOUR, alpha)

    @cached_property
    def direction_totals(self) -> tuple[int, int]:
        """ (towards benign, towards pathogenic) over every reclassification in the window """
        towards_benign = 0
        towards_pathogenic = 0
        for transition in self.significance_transitions:
            moved_benign = ClinicalSignificance.distance(transition.from_clinical_significance,
                                                         transition.to_clinical_significance)
            if moved_benign:
                if moved_benign > 0:
                    towards_benign += transition.count
                else:
                    towards_pathogenic += transition.count
        return towards_benign, towards_pathogenic

    @property
    def towards_benign(self) -> int:
        return self.direction_totals[0]

    @property
    def towards_pathogenic(self) -> int:
        return self.direction_totals[1]

    @property
    def benign_share(self) -> Optional[float]:
        moved = self.towards_benign + self.towards_pathogenic
        return round(100 * self.towards_benign / moved, 1) if moved else None

    @property
    def significance_matrix_totals(self) -> list[int]:
        return [sum(row.cells[index].count for row in self.significance_matrix)
                for index in range(len(SIGNIFICANCE_ORDER))]

    @property
    def significance_matrix_total(self) -> int:
        return sum(self.significance_matrix_totals)

    @staticmethod
    def _translucent(colour: str, alpha: float = 0.45) -> str:
        """ Expands the 4 digit chart colours so the sankey bands can be laid over each other """
        digits = colour.lstrip("#")
        if len(digits) == 3:
            digits = "".join(digit * 2 for digit in digits)
        r, g, b = (int(digits[index:index + 2], 16) for index in (0, 2, 4))
        return f"rgba({r},{g},{b},{alpha})"

    # -- 3. Re-evaluation against reclassification, year by year ---------------------------------

    @cached_property
    def yearly_activity(self) -> list[YearlyActivity]:
        first_event_date = self.events_qs.order_by('reclassified_date') \
            .values_list('reclassified_date', flat=True).first()
        if not first_event_date:
            return []

        reviewed_by_year = self._records_by_year(self.reviews_qs)
        reclassified_by_year = self._records_by_year(self.reclassifications_qs)
        return [
            YearlyActivity(year=year,
                           population=ReclassificationEvent.latest_state_qs(
                               as_of=date(year, 1, 1), base_qs=self.timeline_qs).count(),
                           reviewed=reviewed_by_year.get(year, 0),
                           reclassified=reclassified_by_year.get(year, 0))
            for year in range(first_event_date.year, self.window_end.year + 1)
        ]

    @staticmethod
    def _records_by_year(events_qs: QuerySet[ReclassificationEvent]) -> dict[int, int]:
        """ Records touched per year - a record edited twice in a year is one record that year """
        return {row['year']: row['total'] for row in
                events_qs.annotate(year=ExtractYear('reclassified_date')).values('year')
                .annotate(total=Count('classification_id', distinct=True)).order_by()}

    @cached_property
    def serial_activity(self) -> list[SerialCount]:
        """ Whether the same records keep coming back or the effort spreads across the catalogue """
        reviewed = self._years_touched(self.reviews_qs)
        reclassified = self._years_touched(self.reclassifications_qs)
        rows = [SerialCount(label="Never",
                            reviewed=self.classification_count - sum(reviewed.values()),
                            reclassified=self.classification_count - sum(reclassified.values()))]
        for years in range(1, SERIAL_YEAR_CAP + 1):
            label = f"{years}+ years" if years == SERIAL_YEAR_CAP \
                else ("1 year" if years == 1 else f"{years} years")
            rows.append(SerialCount(label=label, reviewed=reviewed.get(years, 0),
                                    reclassified=reclassified.get(years, 0)))
        return rows

    @staticmethod
    def _years_touched(events_qs: QuerySet[ReclassificationEvent]) -> dict[int, int]:
        touched = defaultdict(int)
        for row in events_qs.values('classification_id') \
                .annotate(years=Count(ExtractYear('reclassified_date'), distinct=True)).order_by():
            touched[min(row['years'], SERIAL_YEAR_CAP)] += 1
        return touched

    # -- 4. Survival curves ----------------------------------------------------------------------

    @cached_property
    def survival_cohort_qs(self) -> QuerySet[ReclassificationEvent]:
        """ Everything sitting at the chosen significance on the origin date, however long it had been there """
        return ReclassificationEvent.latest_state_qs(as_of=self.origin, base_qs=self.timeline_qs) \
            .filter(to_clinical_significance=self.origin_significance)

    @cached_property
    def survival(self) -> Optional[SurvivalAnalysis]:
        cohort_ids = set(self.survival_cohort_qs.values_list('classification_id', flat=True))
        span_days = (self.window_end - self.origin).days
        if not cohort_ids or span_days < DAYS_PER_YEAR:
            return None

        follow_up_qs = self.timeline_qs.filter(
            classification_id__in=self.survival_cohort_qs.values('classification_id'),
            reclassified_date__gt=self.origin, reclassified_date__lte=self.window_end)
        first_review = {row['classification_id']: row['first'] for row in
                        follow_up_qs.filter(event_type__in=ReclassificationEventType.review_types())
                        .values('classification_id').annotate(first=Min('reclassified_date')).order_by()}

        first_reclassification: dict[int, tuple[date, Optional[int]]] = {}
        for classification_id, reclassified_date, delta in follow_up_qs \
                .filter(event_type=ReclassificationEventType.RECLASSIFICATION) \
                .order_by('classification_id', 'reclassified_date') \
                .values_list('classification_id', 'reclassified_date', 'significance_delta'):
            first_reclassification.setdefault(classification_id, (reclassified_date, delta))

        withdrawn_on = self._withdrawal_dates(cohort_ids)
        observed_days = {}
        for classification_id in cohort_ids:
            left_observation = min(withdrawn_on.get(classification_id, self.window_end), self.window_end)
            if left_observation > self.origin:
                # a record already withdrawn on the origin date was never there to be revisited
                observed_days[classification_id] = (left_observation - self.origin).days
        if not observed_days:
            return None

        interval_days = DAYS_PER_YEAR * SURVIVAL_INTERVAL_YEARS
        grid = [round(index * interval_days) for index in range(int(span_days / interval_days) + 1)]

        def members(event_dates: dict[int, date]) -> list[tuple[Optional[int], int]]:
            return [((event_dates[classification_id] - self.origin).days
                     if classification_id in event_dates else None, days)
                    for classification_id, days in observed_days.items()]

        first_review = {classification_id: moved_on for classification_id, moved_on in first_review.items()
                        if classification_id in observed_days}
        first_reclassification = {classification_id: moved
                                  for classification_id, moved in first_reclassification.items()
                                  if classification_id in observed_days}
        reviewed_members = members(first_review)
        reclassified_members = members({classification_id: moved_on for classification_id, (moved_on, _)
                                        in first_reclassification.items()})

        return SurvivalAnalysis(
            cohort_size=len(observed_days),
            origin=self.origin,
            window_end=self.window_end,
            origin_significance=self.origin_significance,
            intervals=[round(index * SURVIVAL_INTERVAL_YEARS, 1) for index in range(len(grid))],
            reviewed=SurvivalCurve(label="Not yet re-evaluated", event_label="Re-evaluated",
                                   colour=REVIEW_COLOUR,
                                   survival=self._life_table(reviewed_members, grid),
                                   events=len(first_review), censored=len(withdrawn_on)),
            reclassified=SurvivalCurve(label="Not yet reclassified", event_label="Reclassified",
                                       colour=CHANGE_COLOUR,
                                       survival=self._life_table(reclassified_members, grid),
                                       events=len(first_reclassification), censored=len(withdrawn_on)))

    @property
    def survival_summary(self) -> list[SurvivalCurve]:
        if survival := self.survival:
            return [survival.reviewed, survival.reclassified]
        return []

    def _withdrawal_dates(self, cohort_ids: set[int]) -> dict[int, date]:
        """ A withdrawn record stops being observable, so it leaves the risk set when its flag went up """
        classification_by_collection = dict(
            Classification.objects.filter(pk__in=cohort_ids, withdrawn=True,
                                          flag_collection__isnull=False)
            .values_list('flag_collection_id', 'pk'))
        if not classification_by_collection:
            return {}
        raised = Flag.objects.filter(collection_id__in=classification_by_collection,
                                     flag_type=classification_flag_types.classification_withdrawn) \
            .values('collection_id').annotate(first=Min('created')).order_by()
        return {classification_by_collection[row['collection_id']]: localdate(row['first'])
                for row in raised}

    @staticmethod
    def _life_table(members: list[tuple[Optional[int], int]], grid: list[int]) -> list[float]:
        """
        Share of the cohort still waiting for the event at each grid point, censoring the members that left
        observation part way through the interval they left in.
        :param members: (days until the event, days until the record left observation) per cohort member
        :param grid: interval end points in days, starting at 0
        """
        events_in = [0] * len(grid)
        censored_in = [0] * len(grid)

        def interval_of(days: int) -> int:
            return min(max(bisect_left(grid, days), 1), len(grid) - 1)

        for event_days, censor_days in members:
            if event_days is not None and event_days <= censor_days:
                events_in[interval_of(event_days)] += 1
            else:
                censored_in[interval_of(censor_days)] += 1

        survival = [1.0]
        at_risk = len(members)
        remaining = 1.0
        for index in range(1, len(grid)):
            events, censored = events_in[index], censored_in[index]
            exposed = at_risk - censored / 2
            if events and exposed > 0:
                remaining *= 1 - events / exposed
            survival.append(round(remaining, 4))
            at_risk -= events + censored
        return survival

    # -- 5. Lab league table ---------------------------------------------------------------------

    @cached_property
    def lab_activity(self) -> list[LabActivity]:
        held = self._records_by_lab(
            ReclassificationEvent.latest_state_qs(as_of=self.window_end, base_qs=self.timeline_qs))
        if not held:
            return []
        reviewed = self._records_by_lab(self.reviews_qs)
        reclassified = self._records_by_lab(self.reclassifications_qs)

        active_rates = [reclassified.get(lab_id, 0) / held[lab_id] for lab_id in held
                        if held[lab_id] and (reviewed.get(lab_id) or reclassified.get(lab_id))]
        median_rate = statistics.median(active_rates) if active_rates else 0.0

        labs = {lab.pk: lab for lab in Lab.objects.filter(pk__in=held).select_related('organization')}
        activity = [
            LabActivity(lab_name=labs[lab_id].name,
                        organization_name=labs[lab_id].organization.shortest_name,
                        held=held[lab_id],
                        reviewed=reviewed.get(lab_id, 0),
                        reclassified=reclassified.get(lab_id, 0),
                        baseline=round(median_rate * held[lab_id], 1),
                        is_selected=bool(self.lab and self.lab.pk == lab_id))
            for lab_id in held if lab_id in labs
        ]
        return sorted(activity, key=lambda row: -row.held)

    @staticmethod
    def _records_by_lab(events_qs: QuerySet[ReclassificationEvent]) -> dict[int, int]:
        return {row['lab']: row['total'] for row in
                events_qs.values('lab').annotate(total=Count('classification_id', distinct=True)).order_by()}

    @property
    def labs_with_no_reclassification(self) -> int:
        return len([row for row in self.lab_activity if not row.reclassified])

    # -- 4. Time to reclassification ------------------------------------------------------------

    @cached_property
    def time_to_reclassification(self) -> list[Distribution]:
        initial_date = ReclassificationEvent.objects \
            .filter(classification_id=OuterRef('classification_id'),
                    event_type=ReclassificationEventType.INITIAL) \
            .values('reclassified_date')[:1]

        days_by_significance = defaultdict(list)
        for from_cs, reclassified_date, initial in self.reclassifications_qs \
                .annotate(initial_date=Subquery(initial_date)) \
                .values_list('from_clinical_significance', 'reclassified_date', 'initial_date'):
            if initial:
                days_by_significance[from_cs].append((reclassified_date - initial).days)

        return [
            Distribution(label=significance_label(cs), clinical_significance=cs, days=days_by_significance[cs])
            for cs in SIGNIFICANCE_ORDER if days_by_significance[cs]
        ]

    @cached_property
    def time_to_reclassification_bins(self) -> dict[str, Any]:
        """
        Counts per log spaced bucket of days+1, so the first few months get the room they deserve.
        Positions are log10 so the template can plot them against a plain linear axis.
        """
        distributions = self.time_to_reclassification
        longest = max((max(distribution.days) for distribution in distributions), default=0)
        if not distributions:
            return {"edges": [], "series": [], "tick_values": [], "tick_labels": []}

        span = log10(longest + 2)
        bin_count = max(1, ceil(span * LOG_BINS_PER_DECADE))
        width = span / bin_count
        edges = [index * width for index in range(bin_count + 1)]

        series = []
        for distribution in distributions:
            counts = [0] * bin_count
            for days in distribution.days:
                counts[min(int(log10(days + 1) / width), bin_count - 1)] += 1
            series.append({"label": f"from {distribution.label}", "colour": distribution.colour,
                           "counts": counts})

        ticks = [(days, label) for days, label in DAY_TICKS if days <= longest] or [DAY_TICKS[0]]
        return {
            "edges": edges,
            "series": series,
            "tick_values": [log10(days + 1) for days, _ in ticks],
            "tick_labels": [label for _, label in ticks],
        }

    # -- 5. VUS burden by gene ------------------------------------------------------------------

    @cached_property
    def gene_burden(self) -> list[GeneBurden]:
        latest_qs = ReclassificationEvent.latest_state_qs(as_of=self.date_to, base_qs=self.timeline_qs) \
            .filter(gene_symbol__isnull=False)

        totals = {row['gene_symbol']: row['total'] for row in
                  latest_qs.values('gene_symbol').annotate(total=Count('pk')).order_by()}
        vus_counts = latest_qs.filter(to_clinical_significance=ClinicalSignificance.VUS) \
            .values('gene_symbol').annotate(total=Count('pk')).order_by('-total')[:TOP_GENE_COUNT]

        return [
            GeneBurden(gene_symbol=row['gene_symbol'], vus_count=row['total'],
                       total_count=totals.get(row['gene_symbol'], row['total']))
            for row in vus_counts
        ]

    # -- 6. Evidence driving reclassification ---------------------------------------------------

    @cached_property
    def evidence_key_order(self) -> list[EvidenceKey]:
        """ Criteria lead in ACMG strength order, then the rest as they appear on a classification """
        evidence_keys = EvidenceKeyMap.cached()
        criteria = sorted(evidence_keys.criteria(),
                          key=lambda e_key: (-CriteriaEvaluation.POINTS.get(e_key.default_crit_evaluation, 0),
                                             e_key.pretty_label.lower()))
        criteria_keys = {e_key.key for e_key in criteria}
        return criteria + [e_key for e_key in evidence_keys.all_keys
                           if e_key.key not in criteria_keys
                           and e_key.key not in EVIDENCE_KEYS_EXCLUDED]

    @cached_property
    def criteria_keys(self) -> set[str]:
        return {e_key.key for e_key in EvidenceKeyMap.cached().criteria()}

    @cached_property
    def evidence_group_by_key(self) -> dict[str, str]:
        """
        Criteria are namespaced per lab config, so a cross lab view folds each lab's version of a criterion
        together while a single lab view keeps that lab's own keys.
        """
        if self.lab:
            return {}
        return {e_key.key: e_key.without_namespace for e_key in EvidenceKeyMap.cached().criteria()}

    def evidence_group(self, key: str) -> str:
        return self.evidence_group_by_key.get(key, key)

    @cached_property
    def evidence_group_order(self) -> list[tuple[str, str]]:
        """ (group, label) in display order, labelled by the first key that lands in each group """
        groups = []
        seen = set()
        for e_key in self.evidence_key_order:
            group = self.evidence_group(e_key.key)
            if group not in seen:
                seen.add(group)
                groups.append((group, e_key.pretty_label))
        return groups

    @cached_property
    def evidence_towards_pathogenic(self) -> list[EvidenceMovement]:
        return self._evidence_movements(towards_benign=False)

    @cached_property
    def evidence_towards_benign(self) -> list[EvidenceMovement]:
        return self._evidence_movements(towards_benign=True)

    def _evidence_movements(self, towards_benign: bool) -> list[EvidenceMovement]:
        applied = defaultdict(int)
        unapplied = defaultdict(int)
        strengthened = defaultdict(int)
        weakened = defaultdict(int)
        changed = defaultdict(int)

        for key, row_towards_benign, old_value, new_value, total in self._evidence_movement_counts:
            if row_towards_benign != towards_benign:
                continue
            group = self.evidence_group(key)
            if key not in self.criteria_keys:
                changed[group] += total
                continue

            old_met = CriteriaEvaluation.is_met(old_value)
            new_met = CriteriaEvaluation.is_met(new_value)
            if old_met != new_met:
                (applied if new_met else unapplied)[group] += total
            else:
                old_strength = abs(CriteriaEvaluation.POINTS.get(old_value, 0))
                new_strength = abs(CriteriaEvaluation.POINTS.get(new_value, 0))
                if new_strength > old_strength:
                    strengthened[group] += total
                elif new_strength < old_strength:
                    weakened[group] += total
                else:
                    changed[group] += total

        moved = set(applied) | set(unapplied) | set(strengthened) | set(weakened) | set(changed)
        criteria_groups = {self.evidence_group(key) for key in self.criteria_keys}
        return [EvidenceMovement(key=group, label=label, applied=applied[group], unapplied=unapplied[group],
                                 strengthened=strengthened[group], weakened=weakened[group],
                                 changed=changed[group], is_criteria=group in criteria_groups)
                for group, label in self.evidence_group_order if group in moved]

    @staticmethod
    def _collapse_evidence(movements: list[EvidenceMovement]) -> tuple[list[EvidenceMovement], int]:
        """ Criteria and the busiest few other keys keep their own row, the long tail folds into one """
        others = sorted((movement for movement in movements if not movement.is_criteria),
                        key=lambda movement: -movement.total)
        kept = {movement.key for movement in others[:TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT]}
        folded = others[TOP_NON_CRITERIA_EVIDENCE_KEY_COUNT:]
        rows = [movement for movement in movements if movement.is_criteria or movement.key in kept]
        if folded:
            rows.append(EvidenceMovement(
                key=OTHER_EVIDENCE_GROUP,
                label=f"Other ({len(folded)} keys)",
                applied=sum(movement.applied for movement in folded),
                unapplied=sum(movement.unapplied for movement in folded),
                strengthened=sum(movement.strengthened for movement in folded),
                weakened=sum(movement.weakened for movement in folded),
                changed=sum(movement.changed for movement in folded),
                is_criteria=False))
        return rows, len(folded)

    @cached_property
    def _evidence_movement_counts(self) -> list[tuple[str, bool, Optional[str], Optional[str], int]]:
        event_ids = list(self.reclassifications_qs.filter(from_modification__isnull=False)
                         .values_list('pk', flat=True))
        if not event_ids:
            return []

        cache_key = f"reclassification_evidence_movements:{self.query_string}:{len(event_ids)}"
        counts = cache.get(cache_key)
        if counts is None:
            counts = self._count_evidence_movements(event_ids, sorted(self.criteria_keys))
            cache.set(cache_key, counts, EVIDENCE_MOVEMENT_CACHE_SECONDS)
        return counts

    @staticmethod
    def _count_evidence_movements(event_ids: list[int], criteria_keys: list[str]) \
            -> list[tuple[str, bool, Optional[str], Optional[str], int]]:
        """
        Diffs each event's two published_evidence blobs in the database - shipping thousands of whole
        evidence documents back to python to compare them would be far more expensive. Criteria carry the
        strength on each side so the caller can tell an on/off from a strength shift, everything else only
        needs to say that it moved.
        """
        event_table = ReclassificationEvent._meta.db_table
        modification_table = ClassificationModification._meta.db_table
        sql = f"""
            SELECT movement.key,
                   movement.towards_benign,
                   movement.old_value,
                   movement.new_value,
                   COUNT(*) AS total
            FROM (
                SELECT event.id,
                       entry.key,
                       event.significance_delta > 0 AS towards_benign,
                       CASE WHEN entry.key = ANY(%s) THEN
                           COALESCE(old_version.published_evidence -> entry.key -> 'value',
                                    old_version.published_evidence -> entry.key) #>> '{{}}' END AS old_value,
                       CASE WHEN entry.key = ANY(%s) THEN
                           COALESCE(new_version.published_evidence -> entry.key -> 'value',
                                    new_version.published_evidence -> entry.key) #>> '{{}}' END AS new_value
                FROM {event_table} event
                JOIN {modification_table} old_version ON old_version.id = event.from_modification_id
                JOIN {modification_table} new_version ON new_version.id = event.to_modification_id
                CROSS JOIN LATERAL (
                    SELECT key FROM jsonb_object_keys(COALESCE(old_version.published_evidence, '{{}}'::jsonb)) AS key
                    UNION
                    SELECT key FROM jsonb_object_keys(COALESCE(new_version.published_evidence, '{{}}'::jsonb)) AS key
                ) entry
                WHERE event.id = ANY(%s)
                  AND event.significance_delta <> 0
                  AND COALESCE(old_version.published_evidence -> entry.key -> 'value',
                               old_version.published_evidence -> entry.key)
                      IS DISTINCT FROM
                      COALESCE(new_version.published_evidence -> entry.key -> 'value',
                               new_version.published_evidence -> entry.key)
            ) movement
            GROUP BY 1, 2, 3, 4
        """
        with connection.cursor() as cursor:
            cursor.execute(sql, [criteria_keys, criteria_keys, event_ids])
            return [(key, towards_benign, old_value, new_value, total)
                    for key, towards_benign, old_value, new_value, total in cursor.fetchall()]

    # -- 7. Points to cross a boundary -----------------------------------------------------------

    @cached_property
    def points_transitions(self) -> list[PointsTransition]:
        """ How far the met ACMG criteria had to travel for a record to cross each bucket boundary """
        deltas = defaultdict(list)
        for from_cs, to_cs, from_points, to_points in self.reclassifications_qs \
                .filter(from_points__isnull=False, to_points__isnull=False) \
                .values_list('from_clinical_significance', 'to_clinical_significance',
                             'from_points', 'to_points'):
            deltas[(from_cs, to_cs)].append(to_points - from_points)

        transitions = [PointsTransition(from_clinical_significance=from_cs, to_clinical_significance=to_cs,
                                        deltas=points)
                       for (from_cs, to_cs), points in deltas.items()]
        return sorted(transitions, key=lambda transition: -transition.count)

    # -- page furniture -------------------------------------------------------------------------

    @property
    def _activity_chart(self) -> dict[str, Any]:
        # a year that opened with nothing in the catalogue has no rate to plot, only a row in the table
        activity = [year for year in self.yearly_activity if year.population]
        return {
            "years": [year.year for year in activity],
            "reviewed_percents": [year.reviewed_percent for year in activity],
            "reclassified_percents": [year.reclassified_percent for year in activity],
            "reviews_per_change": [year.reviews_per_change for year in activity],
            "reviewed": [year.reviewed for year in activity],
            "reclassified": [year.reclassified for year in activity],
            "populations": [year.population for year in activity],
            "serial_labels": [row.label for row in self.serial_activity],
            "serial_reviewed": [row.reviewed for row in self.serial_activity],
            "serial_reclassified": [row.reclassified for row in self.serial_activity],
            "review_colour": REVIEW_COLOUR,
            "change_colour": CHANGE_COLOUR,
        }

    @property
    def _survival_chart(self) -> dict[str, Any]:
        survival = self.survival
        if not survival:
            return {"intervals": []}
        return {
            "intervals": survival.intervals,
            "curves": [{"label": curve.label, "colour": curve.colour, "survival": curve.survival}
                       for curve in (survival.reviewed, survival.reclassified)],
        }

    @property
    def _lab_chart(self) -> dict[str, Any]:
        return {
            "labs": [row.lab_name for row in self.lab_activity],
            "reviewed_percents": [row.reviewed_percent for row in self.lab_activity],
            "reclassified_percents": [row.reclassified_percent for row in self.lab_activity],
            "held": [row.held for row in self.lab_activity],
            "review_colour": REVIEW_COLOUR,
            "change_colour": CHANGE_COLOUR,
        }

    @property
    def _points_chart(self) -> dict[str, Any]:
        return {
            "labels": [transition.label for transition in self.points_transitions],
            "deltas": [transition.deltas for transition in self.points_transitions],
            "colours": [transition.colour for transition in self.points_transitions],
        }

    @cached_property
    def chart_data(self) -> dict[str, Any]:
        """ Everything the plotly calls need, in one JSON blob """
        return {
            "flow": self.significance_flow,
            "activity": self._activity_chart,
            "survival": self._survival_chart,
            "labs": self._lab_chart,
            "points": self._points_chart,
            "time_to_reclassification": self.time_to_reclassification_bins,
            "gene_burden": {
                "symbols": [burden.gene_symbol for burden in self.gene_burden],
                "vus_counts": [burden.vus_count for burden in self.gene_burden],
                "percents": [burden.percent for burden in self.gene_burden],
                "colour": NEUTRAL_COUNT_COLOUR,
            },
            "evidence_towards_pathogenic": self._evidence_chart(self.evidence_towards_pathogenic),
            "evidence_towards_benign": self._evidence_chart(self.evidence_towards_benign),
        }

    def _evidence_chart(self, movements: list[EvidenceMovement]) -> dict[str, Any]:
        """ The collapsed rows are what's drawn, the expanded rows are what the other row opens up to """
        collapsed, folded_count = self._collapse_evidence(movements)
        return {
            "collapsed": self._evidence_series(collapsed),
            "expanded": self._evidence_series(movements),
            "other_key": OTHER_EVIDENCE_GROUP,
            "folded_count": folded_count,
            "applied_colour": EVIDENCE_APPLIED_COLOUR,
            "unapplied_colour": EVIDENCE_UNAPPLIED_COLOUR,
            "strengthened_colour": EVIDENCE_STRENGTHENED_COLOUR,
            "weakened_colour": EVIDENCE_WEAKENED_COLOUR,
            "changed_colour": EVIDENCE_CHANGED_COLOUR,
        }

    @staticmethod
    def _evidence_series(movements: list[EvidenceMovement]) -> dict[str, Any]:
        return {
            "keys": [movement.key for movement in movements],
            "labels": [movement.label for movement in movements],
            "applied": [movement.applied for movement in movements],
            "unapplied": [-movement.unapplied for movement in movements],
            "strengthened": [movement.strengthened for movement in movements],
            "weakened": [-movement.weakened for movement in movements],
            "changed": [movement.changed for movement in movements],
        }

    @property
    def significance_legend(self) -> list[dict[str, str]]:
        return [{"label": significance_label(cs), "css_class": ClinicalSignificance.css_class(cs)}
                for cs in SIGNIFICANCE_ORDER]

    @property
    def query_string(self) -> str:
        parts = []
        if self.lab:
            parts.append(f"lab={self.lab.pk}")
        elif self.organization:
            parts.append(f"organization={self.organization.pk}")
        if self.date_from:
            parts.append(f"date_from={self.date_from.isoformat()}")
        if self.date_to:
            parts.append(f"date_to={self.date_to.isoformat()}")
        parts.append(f"lab_external={self.lab_external.value}")
        parts.append(f"origin={self.origin.isoformat()}")
        parts.append(f"origin_significance={self.origin_significance}")
        return "&".join(parts)


def _timelines_up_to_date() -> ReclassificationBuildResult:
    """ Opening the page is what keeps the timelines current, unless there's more waiting than a
        web worker should take on - then the nightly task's queue picks the batch up. """
    build_result = ReclassificationEventBuilder.bring_up_to_date(
        max_classifications=RECLASSIFICATION_PAGE_BUILD_LIMIT)
    if build_result.outstanding:
        reclassification_events_update.delay()
    return build_result


@require_superuser
def view_reclassification_analytics(request: HttpRequest) -> HttpResponseBase:
    build_result = _timelines_up_to_date()
    analytics = ReclassificationAnalytics.from_request(request)
    return render(request, "classification/classification_reclassification_analytics.html", {
        "analytics": analytics,
        "build_result": build_result,
        "lab_external_choices": LabExternalFilter.choices,
        "lab_external_default": analytics.lab_external,
        "significance_choices": [(cs, significance_label(cs)) for cs in SIGNIFICANCE_ORDER],
    })
