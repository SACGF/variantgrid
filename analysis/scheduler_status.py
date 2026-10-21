"""
The node scheduler's live queue as an operator sees it, for the Analysis Status page
(`analysis/views/views.py:view_analysis_status`).

Entry points: `get_scheduler_status` (deployment-wide counts plus the broker depths of the two
scheduler queues) and `get_analyses_with_work` (the same counts per analysis). Every bucket is a
condition the scheduler itself acts on in `analysis/tasks/analysis_update_tasks.py` - a live lease
is `_live_lease_q`, backing off is the `run_after` gate in `lease_ready_nodes`, and lapsed is what
`lease_ready_nodes` reclaims (QUEUED) or permanently fails (LOADING / LOADING_CACHE) on its next
pass - so the page and the dispatcher cannot disagree.
"""
from dataclasses import dataclass
from datetime import datetime

from django.conf import settings
from django.db.models import Count, Exists, Min, OuterRef, Q, QuerySet, Subquery

from analysis.enums import NodeStatus
from analysis.models.nodes.analysis_node import AnalysisNode
from analysis.models.nodes.node_version import NodeTask
from analysis.tasks.analysis_update_tasks import _live_lease_q
from library.vg.status import broker_queue_depths

SCHEDULER_QUEUES = ["analysis_workers", "scheduling_single_worker"]
WORKER_STATUSES = [NodeStatus.LOADING, NodeStatus.LOADING_CACHE]


@dataclass(frozen=True)
class SchedulerStatus:
    waiting: int  # DIRTY, not backing off
    backing_off: int  # DIRTY, current-version NodeTask has run_after in the future
    dispatched: int  # QUEUED with a live lease - handed to analysis_workers, not yet started
    running: int  # LOADING / LOADING_CACHE with a live lease
    lapsed: int  # QUEUED / LOADING / LOADING_CACHE whose current-version lease is missing or expired
    analyses_with_work: int
    in_flight_target: int
    queues: dict[str, dict]  # broker_queue_depths of SCHEDULER_QUEUES

    @property
    def in_flight(self) -> int:
        return self.dispatched + self.running

    @property
    def graph_counts(self) -> dict[str, int]:
        """ Series name -> count for the page's stacked bar, in stacking order """
        return {"Waiting": self.waiting, "Backing off": self.backing_off, "Dispatched": self.dispatched,
                "Running": self.running, "Lapsed": self.lapsed}

    @property
    def has_work(self) -> bool:
        return any(self.graph_counts.values())


@dataclass(frozen=True)
class AnalysisWork:
    analysis_id: int
    analysis_name: str
    username: str
    waiting: int  # DIRTY, backing off included
    dispatched: int
    running: int
    lapsed: int
    oldest_modified: datetime  # oldest `modified` among the analysis's loading nodes


def _loading_nodes(now) -> QuerySet:
    """ One row per loading node (subqueries, not joins), so conditional counts need no distinct """
    current_task = NodeTask.objects.filter(node_version__node_id=OuterRef("pk"),
                                           node_version__version=OuterRef("version"))
    live_lease = AnalysisNode.objects.filter(pk=OuterRef("pk")).filter(_live_lease_q(now))
    return (AnalysisNode.objects.filter(status__in=NodeStatus.LOADING_STATUSES)
            .annotate(live_lease=Exists(live_lease),
                      current_run_after=Subquery(current_task.values("run_after")[:1])))


def _dirty_counts(now) -> dict[str, Count]:
    not_backing_off = Q(current_run_after__isnull=True) | Q(current_run_after__lte=now)
    return {
        "waiting": Count("pk", filter=Q(status=NodeStatus.DIRTY) & not_backing_off),
        "backing_off": Count("pk", filter=Q(status=NodeStatus.DIRTY, current_run_after__gt=now)),
    }


def _leased_counts() -> dict[str, Count]:
    return {
        "dispatched": Count("pk", filter=Q(status=NodeStatus.QUEUED, live_lease=True)),
        "running": Count("pk", filter=Q(status__in=WORKER_STATUSES, live_lease=True)),
        "lapsed": Count("pk", filter=Q(status__in=[NodeStatus.QUEUED, *WORKER_STATUSES], live_lease=False)),
    }


def get_scheduler_status(now) -> SchedulerStatus:
    counts = _loading_nodes(now).aggregate(**_dirty_counts(now), **_leased_counts(),
                                           analyses_with_work=Count("analysis_id", distinct=True))
    return SchedulerStatus(**counts,
                           in_flight_target=settings.ANALYSIS_NODE_DISPATCH_BACKLOG_IN_FLIGHT_TARGET,
                           queues=broker_queue_depths(SCHEDULER_QUEUES))


def get_analyses_with_work(now, limit=50) -> list[AnalysisWork]:
    """ Oldest first - the analysis that has been waiting longest is the one an operator looks at """
    rows = (_loading_nodes(now)
            .values("analysis_id", "analysis__name", "analysis__user__username")
            .annotate(waiting=Count("pk", filter=Q(status=NodeStatus.DIRTY)), **_leased_counts(),
                      oldest_modified=Min("modified"))
            .order_by("oldest_modified", "analysis_id")[:limit])
    return [AnalysisWork(analysis_id=row["analysis_id"], analysis_name=row["analysis__name"],
                         username=row["analysis__user__username"], waiting=row["waiting"],
                         dispatched=row["dispatched"], running=row["running"], lapsed=row["lapsed"],
                         oldest_modified=row["oldest_modified"])
            for row in rows]
