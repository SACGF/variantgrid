# NOTE it would be ideal for this to react to notification instead of having to be called directly
from datetime import datetime, timedelta, timezone
from typing import Dict, List

from django.db.models import Q

from flags.models import FlagType, FlagComment, Flag, FlagStatus


class FlagDelta:

    def __init__(self, flag_type: FlagType):
        self.flag_type = flag_type
        self.added = 0
        self.resolved = 0

    def __lt__(self, other):
        return self.flag_type.label < other.flag_type.label


def flag_chanced_since(since) -> List[FlagDelta]:
    #min_age = datetime.utcnow().replace(tzinfo=timezone.utc) - timedelta(
    #    minutes=2)  # give records 2 minutes to matching properly before reporting
    time_range_q = Q(created__gte=since)  # & Q(created__lte=min_age)

    new_comments = FlagComment.objects.filter(time_range_q)
    # find all flags that have had comments added in the time range
    modified_flags = Flag.objects.filter(pk__in=new_comments.values_list("flag_id", flat=True))
    # now want to exclude flags that were created and then closed in that time
    modified_flags = modified_flags.exclude(
        Q(created__gte=since) & Q(resolution__status__in=[FlagStatus.CLOSED, FlagStatus.REJECTED]))\
            .select_related("resolution", "flag_type")

    # also need to filter out flags that were open before since and are still open now (likewise for closed)
    # to know that we need to find the most recent flag comment from before closing, and the last flag comment
    # just loop through all the possibile candidates for now, and then see if things start to get too wasteful
    flag_deltas: Dict[FlagType, FlagDelta] = dict()
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

    def is_open_status(flag_status: FlagStatus):
        # consider closed and rejected the same
        return flag_status == FlagStatus.OPEN

    for flag in modified_flags:
        is_flag_open_now = is_open_status(flag.resolution.status)
        # grab latest comment from before since
        if reference_comment := FlagComment.objects.filter(flag=flag, created__lte=since).select_related('resolution').order_by('-created').first():
            if is_open_status(reference_comment.resolution.status) != is_flag_open_now:
                # flag existed before since date, see if the status has changed after the since date
                count_delta(flag.flag_type, 1 if is_flag_open_now else -1)
        else:
            # no reference comment, so flag was created after the since date
            # don't count it if we created and closed it within the window
            if is_flag_open_now:
                count_delta(flag.flag_type, 1)

    return sorted(flag_deltas.values())
