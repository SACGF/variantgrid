import statistics
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from math import ceil, log10
from functools import cached_property
from typing import Any, Optional

from django.core.cache import cache
from django.db import connection
from django.db.models import Count, OuterRef, QuerySet, Subquery
from django.db.models.functions import TruncMonth
from django.http import HttpRequest
from django.http.response import HttpResponseBase
from django.shortcuts import render
from django.utils.timezone import localdate

from classification.enums import ClinicalSignificance, CriteriaEvaluation, LabExternalFilter, SpecialEKeys
from classification.models import (
    ClassificationModification,
    EvidenceKey,
    EvidenceKeyMap,
    ReclassificationBuildResult,
    ReclassificationEvent,
    ReclassificationEventBuilder,
)
from classification.tasks.classification_reclassification_tasks import reclassification_events_update
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

RECLASSIFICATION_PAGE_BUILD_LIMIT = 5000

TOP_GENE_COUNT = 30
EVIDENCE_MOVEMENT_CACHE_SECONDS = 60 * 60
# the significance itself changes by definition, and the dates get stamped alongside it
EVIDENCE_KEYS_ALWAYS_CHANGED = {
    SpecialEKeys.CLINICAL_SIGNIFICANCE,
    SpecialEKeys.CURATION_DATE,
    SpecialEKeys.CURATION_VERIFIED_DATE,
}
# everything else in a criteria dropdown means the criterion is applied at some strength
CRITERIA_UNMET_VALUES = [CriteriaEvaluation.NOT_MET, CriteriaEvaluation.NOT_APPLICABLE,
                         CriteriaEvaluation.NEUTRAL, ""]
EVIDENCE_APPLIED_COLOUR = "#4a8f5b"
EVIDENCE_UNAPPLIED_COLOUR = "#b5564e"
EVIDENCE_CHANGED_COLOUR = "#4f7ca6"

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
class YearlyRate:
    year: int
    reclassified: int
    population: int

    @property
    def percent(self) -> Optional[float]:
        return round(100 * self.reclassified / self.population, 2) if self.population else None


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
    How often one evidence key moved across a set of reclassifications. Criteria coming on or coming
    off are counted separately, everything else - including a criterion only shifting strength - is
    a plain change.
    """
    key: str
    label: str
    applied: int
    unapplied: int
    changed: int


class ReclassificationAnalytics:
    """
    Every chart on /classification/reclassification_analytics reads from here, so they all share the
    organisation / lab / date filters the user picked. @see ReclassificationEvent
    """

    def __init__(self, organization: Optional[Organization], lab: Optional[Lab],
                 date_from: Optional[date], date_to: Optional[date], lab_external: LabExternalFilter):
        self.organization = organization
        self.lab = lab
        self.date_from = date_from
        self.date_to = date_to
        self.lab_external = lab_external

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

        return ReclassificationAnalytics(
            organization=organization,
            lab=lab,
            date_from=get_date('date_from'),
            date_to=get_date('date_to'),
            lab_external=lab_external)

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
        return self.events_qs.filter(from_clinical_significance__isnull=False)

    @cached_property
    def reclassification_count(self) -> int:
        return self.reclassifications_qs.count()

    @cached_property
    def classification_count(self) -> int:
        return self.timeline_qs.values('classification_id').distinct().count()

    @cached_property
    def labs(self) -> QuerySet[Lab]:
        qs = Lab.objects.filter(organization=self.organization) if self.organization else Lab.objects.all()
        return qs.order_by('organization__name', 'name')

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
                    background=self._heat(counts.get((from_cs, to_cs), 0), busiest, from_cs, to_cs))
                    for to_cs in SIGNIFICANCE_ORDER])
            for from_cs in SIGNIFICANCE_ORDER
        ]

    @staticmethod
    def _heat(count: int, busiest: int, from_cs: str, to_cs: str) -> str:
        """ Square rooted so the quiet cells still read, and capped well short of the pill colours """
        if not count or not busiest:
            return ""
        moved_benign = ClinicalSignificance.distance(from_cs, to_cs) or 0
        colour = BENIGN_DIRECTION_COLOUR if moved_benign > 0 else PATHOGENIC_DIRECTION_COLOUR
        alpha = round(0.06 + 0.34 * (count / busiest) ** 0.5, 3)
        return ReclassificationAnalytics._translucent(colour, alpha)

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

    # -- 2. Trend over time ---------------------------------------------------------------------

    @cached_property
    def monthly_direction(self) -> dict[str, Any]:
        benign_by_month = defaultdict(int)
        pathogenic_by_month = defaultdict(int)
        for row in self.reclassifications_qs \
                .annotate(month=TruncMonth('reclassified_date')) \
                .values('month', 'significance_delta') \
                .annotate(total=Count('pk')).order_by('month'):
            delta = row['significance_delta']
            if delta is None or delta == 0:
                continue
            by_month = benign_by_month if delta > 0 else pathogenic_by_month
            by_month[row['month']] += row['total']

        months = sorted(set(benign_by_month) | set(pathogenic_by_month))
        benign = [benign_by_month[month] for month in months]
        pathogenic = [pathogenic_by_month[month] for month in months]

        return {
            "months": [month.isoformat() for month in months],
            "benign": benign,
            "pathogenic": pathogenic,
            "benign_rolling": self._rolling_average(benign),
            "pathogenic_rolling": self._rolling_average(pathogenic),
            "benign_colour": BENIGN_DIRECTION_COLOUR,
            "pathogenic_colour": PATHOGENIC_DIRECTION_COLOUR,
        }

    @staticmethod
    def _rolling_average(values: list[int], window: int = 12) -> list[Optional[float]]:
        averages: list[Optional[float]] = []
        for index in range(len(values)):
            if index + 1 < window:
                averages.append(None)
            else:
                averages.append(round(sum(values[index + 1 - window:index + 1]) / window, 2))
        return averages

    # -- 3. VUS reclassification rate per year --------------------------------------------------

    @cached_property
    def vus_rates(self) -> list[YearlyRate]:
        first_event_date = self.events_qs.order_by('reclassified_date') \
            .values_list('reclassified_date', flat=True).first()
        if not first_event_date:
            return []
        last_year = (self.date_to or localdate()).year

        rates = []
        for year in range(first_event_date.year, last_year + 1):
            population = ReclassificationEvent.latest_state_qs(
                as_of=date(year, 1, 1), base_qs=self.timeline_qs) \
                .filter(to_clinical_significance=ClinicalSignificance.VUS).count()
            reclassified = self.reclassifications_qs.filter(
                from_clinical_significance=ClinicalSignificance.VUS,
                reclassified_date__year=year).count()
            if population or reclassified:
                rates.append(YearlyRate(year=year, reclassified=reclassified, population=population))
        return rates

    # -- 4. Time to reclassification ------------------------------------------------------------

    @cached_property
    def time_to_reclassification(self) -> list[Distribution]:
        initial_date = ReclassificationEvent.objects \
            .filter(classification_id=OuterRef('classification_id'), step=1) \
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
                           and e_key.key not in EVIDENCE_KEYS_ALWAYS_CHANGED]

    @cached_property
    def evidence_towards_pathogenic(self) -> list[EvidenceMovement]:
        return self._evidence_movements(towards_benign=False)

    @cached_property
    def evidence_towards_benign(self) -> list[EvidenceMovement]:
        return self._evidence_movements(towards_benign=True)

    def _evidence_movements(self, towards_benign: bool) -> list[EvidenceMovement]:
        criteria_keys = {e_key.key for e_key in EvidenceKeyMap.cached().criteria()}
        applied = defaultdict(int)
        unapplied = defaultdict(int)
        changed = defaultdict(int)

        for key, row_towards_benign, old_met, new_met, total in self._evidence_movement_counts:
            if row_towards_benign != towards_benign:
                continue
            if key in criteria_keys and old_met != new_met:
                (applied if new_met else unapplied)[key] += total
            else:
                changed[key] += total

        moved = set(applied) | set(unapplied) | set(changed)
        return [EvidenceMovement(key=e_key.key, label=e_key.pretty_label, applied=applied[e_key.key],
                                 unapplied=unapplied[e_key.key], changed=changed[e_key.key])
                for e_key in self.evidence_key_order if e_key.key in moved]

    @cached_property
    def _evidence_movement_counts(self) -> list[tuple[str, bool, bool, bool, int]]:
        event_ids = list(self.reclassifications_qs.filter(from_modification__isnull=False)
                         .values_list('pk', flat=True))
        if not event_ids:
            return []

        cache_key = f"reclassification_evidence_movements:{self.query_string}:{len(event_ids)}"
        counts = cache.get(cache_key)
        if counts is None:
            counts = self._count_evidence_movements(event_ids)
            cache.set(cache_key, counts, EVIDENCE_MOVEMENT_CACHE_SECONDS)
        return counts

    @staticmethod
    def _count_evidence_movements(event_ids: list[int]) -> list[tuple[str, bool, bool, bool, int]]:
        """
        Diffs each event's two published_evidence blobs in the database - shipping thousands of whole
        evidence documents back to python to compare them would be far more expensive. Criteria carry
        whether each side was met so the caller can tell an on/off from a strength shift.
        """
        event_table = ReclassificationEvent._meta.db_table
        modification_table = ClassificationModification._meta.db_table
        sql = f"""
            SELECT movement.key,
                   movement.towards_benign,
                   COALESCE(movement.old_value <> ALL(%s), FALSE) AS old_met,
                   COALESCE(movement.new_value <> ALL(%s), FALSE) AS new_met,
                   COUNT(*) AS total
            FROM (
                SELECT event.id,
                       entry.key,
                       event.significance_delta > 0 AS towards_benign,
                       COALESCE(old_version.published_evidence -> entry.key -> 'value',
                                old_version.published_evidence -> entry.key) #>> '{{}}' AS old_value,
                       COALESCE(new_version.published_evidence -> entry.key -> 'value',
                                new_version.published_evidence -> entry.key) #>> '{{}}' AS new_value
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
            cursor.execute(sql, [CRITERIA_UNMET_VALUES, CRITERIA_UNMET_VALUES, event_ids])
            return [(key, towards_benign, old_met, new_met, total)
                    for key, towards_benign, old_met, new_met, total in cursor.fetchall()]

    # -- page furniture -------------------------------------------------------------------------

    @property
    def _vus_rate_chart(self) -> dict[str, Any]:
        # a year that opened with no VUSes has no rate to plot, only a row in the table
        rates = [rate for rate in self.vus_rates if rate.population]
        return {
            "years": [rate.year for rate in rates],
            "percents": [rate.percent for rate in rates],
            "reclassified": [rate.reclassified for rate in rates],
            "populations": [rate.population for rate in rates],
            "colour": NEUTRAL_COUNT_COLOUR,
        }

    @cached_property
    def chart_data(self) -> dict[str, Any]:
        """ Everything the plotly calls need, in one JSON blob """
        return {
            "flow": self.significance_flow,
            "monthly": self.monthly_direction,
            "vus_rates": self._vus_rate_chart,
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

    @staticmethod
    def _evidence_chart(movements: list[EvidenceMovement]) -> dict[str, Any]:
        return {
            "labels": [movement.label for movement in movements],
            "applied": [movement.applied for movement in movements],
            "unapplied": [-movement.unapplied for movement in movements],
            "changed": [movement.changed for movement in movements],
            "applied_colour": EVIDENCE_APPLIED_COLOUR,
            "unapplied_colour": EVIDENCE_UNAPPLIED_COLOUR,
            "changed_colour": EVIDENCE_CHANGED_COLOUR,
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
    })
