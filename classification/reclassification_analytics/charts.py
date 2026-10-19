"""
The plotly JSON for /classification/reclassification_analytics. chart_data is the entry point: one blob the
page's showReclassificationAnalytics draws every chart from, built off a ReclassificationAnalytics.
"""
from collections import defaultdict
from math import ceil, log10
from typing import Any

from classification.reclassification_analytics.analytics import ReclassificationAnalytics
from classification.reclassification_analytics.evidence_movement import (
    OTHER_EVIDENCE_GROUP,
    collapse_evidence,
)
from classification.reclassification_analytics.palette import (
    BENIGN_DIRECTION_COLOUR,
    CHANGE_COLOUR,
    EVIDENCE_APPLIED_COLOUR,
    EVIDENCE_CHANGED_COLOUR,
    EVIDENCE_STRENGTHENED_COLOUR,
    EVIDENCE_UNAPPLIED_COLOUR,
    EVIDENCE_WEAKENED_COLOUR,
    NEUTRAL_COUNT_COLOUR,
    PATHOGENIC_DIRECTION_COLOUR,
    REVIEW_COLOUR,
    translucent,
)
from classification.reclassification_analytics.results import (
    SIGNIFICANCE_ORDER,
    EvidenceMovement,
    significance_label,
)
from snpdb.models.models_enums import ClinicalSignificance

# the days axis is logarithmic - most records move within a few months, and a year either side of the
# two year mark is not a distinction anybody reads off a chart
LOG_BINS_PER_DECADE = 8
DAY_TICKS = [(0, "same day"), (7, "1w"), (30, "1m"), (90, "3m"), (182, "6m"), (365, "1y"),
             (730, "2y"), (1826, "5y"), (3652, "10y")]


def chart_data(analytics: ReclassificationAnalytics) -> dict[str, Any]:
    """ Everything the plotly calls need, in one JSON blob """
    return {
        "flow": significance_flow(analytics),
        "activity": _activity_chart(analytics),
        "survival": _survival_chart(analytics),
        "labs": _lab_chart(analytics),
        "points": _points_chart(analytics),
        "time_to_reclassification": time_to_reclassification_bins(analytics),
        "gene_burden": {
            "symbols": [burden.gene_symbol for burden in analytics.gene_burden],
            "vus_counts": [burden.vus_count for burden in analytics.gene_burden],
            "percents": [burden.percent for burden in analytics.gene_burden],
            "colour": NEUTRAL_COUNT_COLOUR,
        },
        "evidence_towards_pathogenic": _evidence_chart(analytics.evidence.towards_pathogenic),
        "evidence_towards_benign": _evidence_chart(analytics.evidence.towards_benign),
    }


def significance_flow(analytics: ReclassificationAnalytics) -> dict[str, Any]:
    """ The sankey """
    # a source column then a target column, so a significance that both loses and gains records
    # shows up on each side rather than looping back on itself
    source_index = {cs: index for index, cs in enumerate(SIGNIFICANCE_ORDER)}
    target_index = {cs: index + len(SIGNIFICANCE_ORDER) for index, cs in enumerate(SIGNIFICANCE_ORDER)}

    sources, targets, values, link_colours = [], [], [], []
    source_totals = defaultdict(int)
    target_totals = defaultdict(int)
    for transition in analytics.significance_transitions:
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
        link_colours.append(translucent(
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


def time_to_reclassification_bins(analytics: ReclassificationAnalytics) -> dict[str, Any]:
    """
    Counts per log spaced bucket of days+1, so the first few months get the room they deserve.
    Positions are log10 so the template can plot them against a plain linear axis.
    """
    distributions = analytics.time_to_reclassification
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


def _activity_chart(analytics: ReclassificationAnalytics) -> dict[str, Any]:
    # a year that opened with nothing in the catalogue has no rate to plot, only a row in the table
    activity = [year for year in analytics.yearly_activity if year.population]
    return {
        "years": [year.year for year in activity],
        "reviewed_percents": [year.reviewed_percent for year in activity],
        "reclassified_percents": [year.reclassified_percent for year in activity],
        "reviews_per_change": [year.reviews_per_change for year in activity],
        "reviewed": [year.reviewed for year in activity],
        "reclassified": [year.reclassified for year in activity],
        "populations": [year.population for year in activity],
        "serial_labels": [row.label for row in analytics.serial_activity],
        "serial_reviewed": [row.reviewed for row in analytics.serial_activity],
        "serial_reclassified": [row.reclassified for row in analytics.serial_activity],
        "review_colour": REVIEW_COLOUR,
        "change_colour": CHANGE_COLOUR,
    }


def _survival_chart(analytics: ReclassificationAnalytics) -> dict[str, Any]:
    survival = analytics.survival
    if not survival:
        return {"intervals": []}
    return {
        "intervals": survival.intervals,
        "curves": [{"label": curve.label, "colour": curve.colour, "survival": curve.survival}
                   for curve in (survival.reviewed, survival.reclassified)],
    }


def _lab_chart(analytics: ReclassificationAnalytics) -> dict[str, Any]:
    return {
        "labs": [row.lab_name for row in analytics.lab_activity],
        "reviewed_percents": [row.reviewed_percent for row in analytics.lab_activity],
        "reclassified_percents": [row.reclassified_percent for row in analytics.lab_activity],
        "held": [row.held for row in analytics.lab_activity],
        "review_colour": REVIEW_COLOUR,
        "change_colour": CHANGE_COLOUR,
    }


def _points_chart(analytics: ReclassificationAnalytics) -> dict[str, Any]:
    return {
        "labels": [transition.label for transition in analytics.points_transitions],
        "deltas": [transition.deltas for transition in analytics.points_transitions],
        "colours": [transition.colour for transition in analytics.points_transitions],
    }


def _evidence_chart(movements: list[EvidenceMovement]) -> dict[str, Any]:
    """ The collapsed rows are what's drawn, the expanded rows are what the other row opens up to """
    collapsed, folded_count = collapse_evidence(movements)
    return {
        "collapsed": _evidence_series(collapsed),
        "expanded": _evidence_series(movements),
        "other_key": OTHER_EVIDENCE_GROUP,
        "folded_count": folded_count,
        "applied_colour": EVIDENCE_APPLIED_COLOUR,
        "unapplied_colour": EVIDENCE_UNAPPLIED_COLOUR,
        "strengthened_colour": EVIDENCE_STRENGTHENED_COLOUR,
        "weakened_colour": EVIDENCE_WEAKENED_COLOUR,
        "changed_colour": EVIDENCE_CHANGED_COLOUR,
    }


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
