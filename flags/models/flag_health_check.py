# NOTE it would be ideal for this to react to notification instead of having to be called directly
from datetime import datetime
from typing import Dict, List, Optional, Iterable

from django.db.models import Q

from flags.models import FlagType, FlagComment, Flag, FlagStatus, FlagResolution


class FlagDelta:

    def __init__(self, flag_type: FlagType):
        self.flag_type = flag_type
        self.added = 0
        self.resolved = 0

    def __lt__(self, other):
        return self.flag_type.label < other.flag_type.label


def flag_chanced_since(since: datetime, flag_types: Optional[Iterable[FlagType]] = None) -> List[FlagDelta]:
    # TODO take a NOW date
    time_range_q = Q(created__gte=since)  # & Q(created__lte=min_age)
    if not flag_types:
        flag_types = FlagType.objects.all()

    new_comments = FlagComment.objects.filter(time_range_q)
    # find all flags that have had comments added in the time range
    modified_flags = Flag.objects.filter(pk__in=new_comments.values_list("flag_id", flat=True), flag_type__in=flag_types)
    # now want to exclude flags that were created and then closed in that time
    modified_flags = modified_flags.exclude(
        Q(created__gte=since) & Q(resolution__status__in=[FlagStatus.CLOSED, FlagStatus.REJECTED])
    ).select_related("resolution", "flag_type")

    # also need to filter out flags that were open before since and are still open now (likewise for closed)
    # to know that we need to find the most recent flag comment from before closing, and the last flag comment
    # just loop through all the possibile candidates for now, and then see if things start to get too wasteful
    flag_deltas: Dict[FlagType, FlagDelta] = {}

    def count_delta(flag_type: FlagType, diff: int):
        delta = flag_deltas.get(flag_type)
        if not delta:
            delta = FlagDelta(flag_type)
            flag_deltas[flag_type] = delta
        if diff == 1:
            delta.added += 1
        elif diff == -1:
            delta.resolved += 1
        else:
            raise ValueError("Diff must be 1 or -1")

    def is_open_resolution(resolution: FlagResolution):
        # consider closed and rejected the same
        # not just when the flag changed
        return resolution.status == FlagStatus.OPEN

    for flag in modified_flags:
        is_flag_open_now = is_open_resolution(flag.resolution)
        # grab latest comment from before since
        if state_prior_to_since := FlagComment.objects.filter(flag=flag, created__lte=since, resolution__isnull=False).select_related('resolution').order_by('-created').first():
            # TODO It would be better if FlagComment resolution was not nullable and always reflected the state of the flag at the time
            # but instead we get the first comment with a resolution (FlagComments should always have a resolution when they open)
            was_open_then = is_open_resolution(state_prior_to_since.resolution)

            if was_open_then != is_flag_open_now:
                # flag existed before since date, see if the status has changed after the since date
                count_delta(flag.flag_type, 1 if is_flag_open_now else -1)
        else:
            # flag didn't exist before since date, so report it if it's open and ignore it if it closed itself
            if is_flag_open_now:
                # flag is currently open and was created within the time window, report it
                count_delta(flag.flag_type, 1)

    return sorted(flag_deltas.values())
