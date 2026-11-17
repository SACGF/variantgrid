"""
Celery pool introspection for the annotation dispatcher (#2667).

Each tool lane's dispatch cap must track the live ``--concurrency`` of the pool it runs on
(``annotation_workers``, or a pipeline's own ``PipelineDef.worker_queue``), set in the systemd units on
servers / the dev ``run_celery`` script, so it can never drift below the real worker count and leave
workers idle. So we derive it from the running pool rather than re-declaring it in settings;
the only setting is ``ANNOTATION_WORKER_SLOTS_FALLBACK`` for when inspection sees no workers.
"""
import logging

from django.conf import settings
from django.core.cache import cache

from library.log_utils import log_traceback
from variantgrid.celery import app

ANNOTATION_WORKERS_QUEUE = "annotation_workers"
# Pool size only changes on worker restart/deploy, and inspect() is a broadcast RPC we don't want on
# every dispatch - so cache it briefly.
_WORKER_SLOTS_CACHE_TTL = 60


def annotation_worker_slots() -> int:
    """ Live size of the annotation_workers pool - @see worker_queue_slots """
    return worker_queue_slots(ANNOTATION_WORKERS_QUEUE)


def worker_queue_slots(queue: str) -> int:
    """ Live size of the pool consuming `queue` (summed --concurrency), cached ~60s. When celery inspection
        sees no workers on it: settings.ANNOTATION_WORKER_SLOTS_FALLBACK for annotation_workers, 1 for any
        other queue (a tool lane on a small pool, eg AnnotSV on heavy_workers, mustn't be oversubscribed). """
    cache_key = f"worker_queue_slots:{queue}"
    slots = cache.get(cache_key)
    if slots is None:
        if slots := live_worker_queue_slots(queue):
            logging.info("worker_queue_slots(%s): live pool = %d", queue, slots)
        else:
            slots = settings.ANNOTATION_WORKER_SLOTS_FALLBACK if queue == ANNOTATION_WORKERS_QUEUE else 1
            logging.info("worker_queue_slots(%s): no workers seen, using fallback %d", queue, slots)
        cache.set(cache_key, slots, _WORKER_SLOTS_CACHE_TTL)
    return slots


def live_worker_queue_slots(queue: str) -> int:
    """ Summed --concurrency of the running workers consuming `queue`: 0 when none answer (or inspection fails) """
    slots = 0
    try:
        inspect = app.control.inspect()
        active_queues = inspect.active_queues() or {}
        stats = inspect.stats() or {}
        for worker_name, queues in active_queues.items():
            if not any(q.get("name") == queue for q in queues):
                continue
            pool = stats.get(worker_name, {}).get("pool", {})
            if concurrency := pool.get("max-concurrency"):
                slots += concurrency
    except Exception:
        log_traceback()
    return slots
