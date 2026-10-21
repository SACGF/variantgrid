"""
The value objects the reclassification analytics page renders - a row of a table, a curve, a distribution -
and the significance order and short labels every chart and table shares. Built by
analytics.py:ReclassificationAnalytics, read by the template and charts.py.
"""
import statistics
from dataclasses import dataclass
from datetime import date
from typing import Optional

from classification.reclassification_analytics.palette import (
    BENIGN_DIRECTION_COLOUR,
    PATHOGENIC_DIRECTION_COLOUR,
)
from classification.reclassification_analytics.survival import SURVIVAL_INTERVAL_YEARS
from snpdb.enums import ClinicalSignificance

# top to bottom / left to right on every chart and table
SIGNIFICANCE_ORDER = [
    ClinicalSignificance.PATHOGENIC,
    ClinicalSignificance.LIKELY_PATHOGENIC,
    ClinicalSignificance.VUS,
    ClinicalSignificance.LIKELY_BENIGN,
    ClinicalSignificance.BENIGN,
]


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
