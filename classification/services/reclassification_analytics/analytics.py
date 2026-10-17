"""
ReclassificationAnalytics: every table on /classification/reclassification_analytics, worked out over the
ReclassificationEvent timelines the organisation / lab / date filters select. Build one with
ReclassificationAnalytics.from_params (the page's GET parameters) or the constructor; charts.py:chart_data
turns one into the page's plotly JSON.
"""
import statistics
from collections import defaultdict
from collections.abc import Mapping
from datetime import date
from functools import cached_property
from typing import Optional

from dateutil.relativedelta import relativedelta
from django.conf import settings
from django.db.models import Count, Min, OuterRef, QuerySet, Subquery
from django.db.models.functions import ExtractYear
from django.utils.timezone import localdate

from classification.enums.classification_enums import (
    ClinicalSignificance,
    LabExternalFilter,
    ReclassificationEventType,
)
from classification.models.classification_reclassification_models import ReclassificationEvent
from classification.services.reclassification_analytics.evidence_movement import EvidenceMovements
from classification.services.reclassification_analytics.palette import (
    CHANGE_COLOUR,
    MATRIX_HEAT_COLOUR,
    REVIEW_COLOUR,
    translucent,
)
from classification.services.reclassification_analytics.results import (
    SIGNIFICANCE_ORDER,
    Distribution,
    GeneBurden,
    LabActivity,
    PointsTransition,
    SerialCount,
    SignificanceMatrixCell,
    SignificanceMatrixRow,
    SignificanceTransition,
    SurvivalAnalysis,
    SurvivalCurve,
    YearlyActivity,
    significance_label,
)
from classification.services.reclassification_analytics.survival import (
    DAYS_PER_YEAR,
    SURVIVAL_INTERVAL_YEARS,
    life_table,
    withdrawal_dates,
)
from snpdb.models.models import Lab, Organization

TOP_GENE_COUNT = 30
SURVIVAL_ORIGIN_YEARS_BACK = 5
SERIAL_YEAR_CAP = 4


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
    def from_params(params: Mapping[str, str]) -> 'ReclassificationAnalytics':
        """ The page's filters, ignoring any value that doesn't parse """
        def get_date(param: str) -> Optional[date]:
            if raw := params.get(param):
                try:
                    return date.fromisoformat(raw)
                except ValueError:
                    pass
            return None

        organization = None
        lab = None
        if organization_id := params.get('organization'):
            organization = Organization.objects.filter(pk=organization_id).first()
        if lab_id := params.get('lab'):
            lab = Lab.objects.filter(pk=lab_id).first()
            if lab:
                organization = lab.organization

        lab_external = LabExternalFilter.ALL
        if settings.CLASSIFICATION_RECLASSIFICATION_ANALYTICS_INTERNAL_ONLY:
            lab_external = LabExternalFilter.INTERNAL
        elif raw_lab_external := params.get('lab_external'):
            if raw_lab_external in LabExternalFilter.values:
                lab_external = LabExternalFilter(raw_lab_external)

        origin_significance = ClinicalSignificance.VUS
        if raw_significance := params.get('origin_significance'):
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
        if self.lab_external == LabExternalFilter.INTERNAL:
            qs = qs.filter(external=False)
        return qs.order_by('organization__name', 'name')

    @cached_property
    def selected_labs(self) -> list[Lab]:
        """ The labs the organisation / lab / internal filters leave holding a timeline """
        lab_ids = self.timeline_qs.values_list('lab', flat=True).distinct()
        return list(Lab.objects.filter(pk__in=lab_ids).select_related('organization')
                    .order_by('organization__name', 'name'))

    @cached_property
    def organizations(self) -> QuerySet[Organization]:
        qs = Organization.objects.filter(active=True)
        if self.lab_external == LabExternalFilter.INTERNAL:
            qs = qs.filter(lab__external=False).distinct()
        return qs.order_by('name')

    # -- 1. Significance transitions --------------------------------------------------------------

    @cached_property
    def significance_transitions(self) -> list[SignificanceTransition]:
        counts = self.reclassifications_qs \
            .values('from_clinical_significance', 'to_clinical_significance') \
            .annotate(total=Count('pk')).order_by('-total')
        return [SignificanceTransition(from_clinical_significance=row['from_clinical_significance'],
                                       to_clinical_significance=row['to_clinical_significance'],
                                       count=row['total']) for row in counts]

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
        return translucent(MATRIX_HEAT_COLOUR, alpha)

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

        withdrawn_on = withdrawal_dates(cohort_ids)
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
                                   survival=life_table(reviewed_members, grid),
                                   events=len(first_review), censored=len(withdrawn_on)),
            reclassified=SurvivalCurve(label="Not yet reclassified", event_label="Reclassified",
                                       colour=CHANGE_COLOUR,
                                       survival=life_table(reclassified_members, grid),
                                       events=len(first_reclassification), censored=len(withdrawn_on)))

    @property
    def survival_summary(self) -> list[SurvivalCurve]:
        if survival := self.survival:
            return [survival.reviewed, survival.reclassified]
        return []

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
    def evidence(self) -> EvidenceMovements:
        return EvidenceMovements(self.reclassifications_qs, fold_lab_namespaces=not self.lab,
                                 cache_scope=self.query_string)

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
