"""
The life table behind the "Waiting for a Second Look" survival curves: life_table does the actuarial maths over
(event, censor) day pairs, withdrawal_dates finds when a withdrawn record left observation.
analytics.py:ReclassificationAnalytics.survival builds the cohort and calls both.
"""
from bisect import bisect_left
from datetime import date
from typing import Optional

from django.db.models import Min
from django.utils.timezone import localdate

from classification.models.classification import Classification
from classification.models.flag_types import classification_flag_types
from flags.models.models import Flag

DAYS_PER_YEAR = 365.25
# most records that move do it inside the first couple of years, so the life table steps every 6 months
SURVIVAL_INTERVAL_YEARS = 0.5


def withdrawal_dates(cohort_ids: set[int]) -> dict[int, date]:
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


def life_table(members: list[tuple[Optional[int], int]], grid: list[int]) -> list[float]:
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
