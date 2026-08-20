import statistics
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
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

from classification.enums import ClinicalSignificance
from classification.models import ClassificationModification, EvidenceKeyMap, ReclassificationEvent
from library.django_utils import require_superuser
from snpdb.models import Lab, Organization

SIGNIFICANCE_ORDER = [
    ClinicalSignificance.BENIGN,
    ClinicalSignificance.LIKELY_BENIGN,
    ClinicalSignificance.VUS,
    ClinicalSignificance.LIKELY_PATHOGENIC,
    ClinicalSignificance.PATHOGENIC,
]

# Diverging benign (blue) - VUS (neutral) - pathogenic (red), matching the .cs- cell colours in global.scss
SIGNIFICANCE_COLOURS = {
    ClinicalSignificance.BENIGN: "#1c5cab",
    ClinicalSignificance.LIKELY_BENIGN: "#5598e7",
    ClinicalSignificance.VUS: "#6f6e6a",
    ClinicalSignificance.LIKELY_PATHOGENIC: "#e8735c",
    ClinicalSignificance.PATHOGENIC: "#b02b2b",
    ClinicalSignificance.OTHER: "#b8b7b2",
}

BENIGN_DIRECTION_COLOUR = "#2a78d6"
PATHOGENIC_DIRECTION_COLOUR = "#b02b2b"
# counts that carry no benign/pathogenic meaning take a hue outside the significance language
NEUTRAL_COUNT_COLOUR = "#4a3aa7"

TOP_GENE_COUNT = 30
TOP_EVIDENCE_KEY_COUNT = 25
EVIDENCE_KEY_DIFF_CACHE_SECONDS = 60 * 60


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
    def colour(self) -> str:
        return SIGNIFICANCE_COLOURS.get(self.clinical_significance)


@dataclass(frozen=True)
class SignificanceTransition:
    from_clinical_significance: str
    to_clinical_significance: str
    count: int

    @property
    def from_label(self) -> str:
        return significance_label(self.from_clinical_significance)

    @property
    def to_label(self) -> str:
        return significance_label(self.to_clinical_significance)

    @property
    def from_colour(self) -> str:
        return SIGNIFICANCE_COLOURS.get(self.from_clinical_significance)

    @property
    def to_colour(self) -> str:
        return SIGNIFICANCE_COLOURS.get(self.to_clinical_significance)


@dataclass(frozen=True)
class YearlyDirection:
    year: int
    benign: int
    pathogenic: int


@dataclass(frozen=True)
class YearlyRate:
    year: int
    reclassified: int
    population: int

    @property
    def percent(self) -> Optional[float]:
        return round(100 * self.reclassified / self.population, 2) if self.population else None

    @property
    def percent_display(self) -> str:
        return f"{self.percent}%" if self.percent is not None else "-"


@dataclass(frozen=True)
class GeneBurden:
    gene_symbol: str
    vus_count: int
    total_count: int

    @property
    def percent(self) -> float:
        return round(100 * self.vus_count / self.total_count, 1) if self.total_count else 0.0


@dataclass(frozen=True)
class EvidenceKeyChange:
    key: str
    label: str
    count: int


class ReclassificationAnalytics:
    """
    Every chart on /classification/reclassification_analytics reads from here, so they all share the
    organisation / lab / date / provenance filters the user picked. @see ReclassificationEvent
    """

    PROVENANCE_LOCAL = "local"
    PROVENANCE_SYNCED = "synced"
    PROVENANCE_ALL = "all"
    PROVENANCE_CHOICES = [
        (PROVENANCE_LOCAL, "Local labs"),
        (PROVENANCE_SYNCED, "Synced labs"),
        (PROVENANCE_ALL, "All labs"),
    ]

    def __init__(self, organization: Optional[Organization], lab: Optional[Lab],
                 date_from: Optional[date], date_to: Optional[date], provenance: str):
        self.organization = organization
        self.lab = lab
        self.date_from = date_from
        self.date_to = date_to
        self.provenance = provenance

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

        provenance = request.GET.get('provenance') or ReclassificationAnalytics.PROVENANCE_LOCAL
        if provenance not in dict(ReclassificationAnalytics.PROVENANCE_CHOICES):
            provenance = ReclassificationAnalytics.PROVENANCE_LOCAL

        return ReclassificationAnalytics(
            organization=organization,
            lab=lab,
            date_from=get_date('date_from'),
            date_to=get_date('date_to'),
            provenance=provenance)

    # -- querysets every chart builds on --------------------------------------------------------

    @cached_property
    def timeline_qs(self) -> QuerySet[ReclassificationEvent]:
        """ The whole timeline for the selected labs, unbounded by date - denominators need the history """
        qs = ReclassificationEvent.objects.all()
        if self.lab:
            qs = qs.filter(lab=self.lab)
        elif self.organization:
            qs = qs.filter(lab__organization=self.organization)
        if self.provenance == ReclassificationAnalytics.PROVENANCE_LOCAL:
            qs = qs.filter(lab__external=False)
        elif self.provenance == ReclassificationAnalytics.PROVENANCE_SYNCED:
            qs = qs.filter(lab__external=True)
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
        for transition in self.significance_transitions:
            from_cs = transition.from_clinical_significance
            to_cs = transition.to_clinical_significance
            if from_cs not in source_index or to_cs not in target_index:
                continue
            sources.append(source_index[from_cs])
            targets.append(target_index[to_cs])
            values.append(transition.count)
            moved_benign = ClinicalSignificance.distance(from_cs, to_cs) or 0
            link_colours.append(self._translucent(
                BENIGN_DIRECTION_COLOUR if moved_benign > 0 else PATHOGENIC_DIRECTION_COLOUR))

        return {
            "labels": [f"{significance_label(cs)} from" for cs in SIGNIFICANCE_ORDER] +
                      [f"{significance_label(cs)} to" for cs in SIGNIFICANCE_ORDER],
            "node_colours": [SIGNIFICANCE_COLOURS[cs] for cs in SIGNIFICANCE_ORDER] * 2,
            "sources": sources,
            "targets": targets,
            "values": values,
            "link_colours": link_colours,
        }

    @staticmethod
    def _translucent(colour: str) -> str:
        r, g, b = (int(colour[index:index + 2], 16) for index in (1, 3, 5))
        return f"rgba({r},{g},{b},0.45)"

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

    @cached_property
    def yearly_direction(self) -> list[YearlyDirection]:
        benign_by_year = defaultdict(int)
        pathogenic_by_year = defaultdict(int)
        for row in self.reclassifications_qs \
                .values('reclassified_date__year', 'significance_delta') \
                .annotate(total=Count('pk')).order_by():
            delta = row['significance_delta']
            if not delta:
                continue
            by_year = benign_by_year if delta > 0 else pathogenic_by_year
            by_year[row['reclassified_date__year']] += row['total']

        return [YearlyDirection(year=year, benign=benign_by_year[year], pathogenic=pathogenic_by_year[year])
                for year in sorted(set(benign_by_year) | set(pathogenic_by_year))]

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

    # -- 6. Evidence keys driving reclassification ----------------------------------------------

    @cached_property
    def evidence_key_changes(self) -> list[EvidenceKeyChange]:
        event_ids = list(self.reclassifications_qs.filter(from_modification__isnull=False)
                         .values_list('pk', flat=True))
        if not event_ids:
            return []

        cache_key = f"reclassification_evidence_key_changes:{self.filter_query}:{len(event_ids)}"
        counts = cache.get(cache_key)
        if counts is None:
            counts = self._count_changed_evidence_keys(event_ids)
            cache.set(cache_key, counts, EVIDENCE_KEY_DIFF_CACHE_SECONDS)

        evidence_keys = EvidenceKeyMap.cached()
        return [
            EvidenceKeyChange(key=key, label=evidence_keys.get(key).pretty_label, count=count)
            for key, count in counts[:TOP_EVIDENCE_KEY_COUNT]
        ]

    @staticmethod
    def _count_changed_evidence_keys(event_ids: list[int]) -> list[tuple[str, int]]:
        """
        Diffs each event's two published_evidence blobs in the database - shipping thousands of whole
        evidence documents back to python to compare them would be far more expensive.
        """
        event_table = ReclassificationEvent._meta.db_table
        modification_table = ClassificationModification._meta.db_table
        sql = f"""
            SELECT changed.key, COUNT(*) AS total
            FROM (
                SELECT event.id, entry.key
                FROM {event_table} event
                JOIN {modification_table} old_version ON old_version.id = event.from_modification_id
                JOIN {modification_table} new_version ON new_version.id = event.to_modification_id
                CROSS JOIN LATERAL (
                    SELECT key FROM jsonb_object_keys(COALESCE(old_version.published_evidence, '{{}}'::jsonb)) AS key
                    UNION
                    SELECT key FROM jsonb_object_keys(COALESCE(new_version.published_evidence, '{{}}'::jsonb)) AS key
                ) entry
                WHERE event.id = ANY(%s)
                  AND COALESCE(old_version.published_evidence -> entry.key -> 'value',
                               old_version.published_evidence -> entry.key)
                      IS DISTINCT FROM
                      COALESCE(new_version.published_evidence -> entry.key -> 'value',
                               new_version.published_evidence -> entry.key)
            ) changed
            GROUP BY changed.key
            ORDER BY total DESC
        """
        with connection.cursor() as cursor:
            cursor.execute(sql, [event_ids])
            return [(key, total) for key, total in cursor.fetchall()]

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
            "colour": SIGNIFICANCE_COLOURS[ClinicalSignificance.VUS],
        }

    @cached_property
    def chart_data(self) -> dict[str, Any]:
        """ Everything the plotly calls need, in one JSON blob """
        return {
            "flow": self.significance_flow,
            "monthly": self.monthly_direction,
            "vus_rates": self._vus_rate_chart,
            "time_to_reclassification": [
                {"label": f"from {distribution.label}", "colour": distribution.colour, "days": distribution.days}
                for distribution in self.time_to_reclassification
            ],
            "gene_burden": {
                "symbols": [burden.gene_symbol for burden in self.gene_burden],
                "vus_counts": [burden.vus_count for burden in self.gene_burden],
                "percents": [burden.percent for burden in self.gene_burden],
                "colour": SIGNIFICANCE_COLOURS[ClinicalSignificance.VUS],
            },
            "evidence_keys": {
                "labels": [change.label for change in self.evidence_key_changes],
                "counts": [change.count for change in self.evidence_key_changes],
                "colour": NEUTRAL_COUNT_COLOUR,
            },
        }


    @property
    def significance_legend(self) -> list[dict[str, str]]:
        return [{"label": significance_label(cs), "colour": SIGNIFICANCE_COLOURS[cs]}
                for cs in SIGNIFICANCE_ORDER]

    @property
    def filter_query(self) -> str:
        parts = []
        if self.lab:
            parts.append(f"lab={self.lab.pk}")
        elif self.organization:
            parts.append(f"organization={self.organization.pk}")
        if self.date_from:
            parts.append(f"date_from={self.date_from.isoformat()}")
        if self.date_to:
            parts.append(f"date_to={self.date_to.isoformat()}")
        parts.append(f"provenance={self.provenance}")
        return "&".join(parts)


@require_superuser
def view_reclassification_analytics(request: HttpRequest) -> HttpResponseBase:
    analytics = ReclassificationAnalytics.from_request(request)
    return render(request, "classification/classification_reclassification_analytics.html", {
        "analytics": analytics,
        "provenance_choices": ReclassificationAnalytics.PROVENANCE_CHOICES,
    })
