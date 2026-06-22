"""
Celery pool introspection for the annotation dispatcher (#2667).

The dispatch cap must track the live ``annotation_workers`` ``--concurrency`` (set in the systemd
units on servers / the dev ``run_celery`` script) so it can never drift below the real worker count
and leave workers idle. So we derive it from the running pool rather than re-declaring it in settings;
the only setting is ``ANNOTATION_WORKER_SLOTS_FALLBACK`` for when inspection sees no workers.
"""
import logging

from django.conf import settings
from django.core.cache import cache

from library.log_utils import log_traceback
from variantgrid.celery import app

ANNOTATION_WORKERS_QUEUE = "annotation_workers"
_WORKER_SLOTS_CACHE_KEY = "annotation_worker_slots"
# Pool size only changes on worker restart/deploy, and inspect() is a broadcast RPC we don't want on
# every dispatch - so cache it briefly.
_WORKER_SLOTS_CACHE_TTL = 60


def annotation_worker_slots() -> int:
    """ Live size of the annotation_workers pool (summed --concurrency), cached ~60s. Falls back to
        settings.ANNOTATION_WORKER_SLOTS_FALLBACK when celery inspection sees no workers. """
    slots = cache.get(_WORKER_SLOTS_CACHE_KEY)
    if slots is None:
        slots = _inspect_annotation_worker_slots()
        cache.set(_WORKER_SLOTS_CACHE_KEY, slots, _WORKER_SLOTS_CACHE_TTL)
    return slots


def _inspect_annotation_worker_slots() -> int:
    fallback = settings.ANNOTATION_WORKER_SLOTS_FALLBACK
    try:
        inspect = app.control.inspect()
        active_queues = inspect.active_queues() or {}
        stats = inspect.stats() or {}
        slots = 0
        for worker_name, queues in active_queues.items():
            if not any(q.get("name") == ANNOTATION_WORKERS_QUEUE for q in queues):
                continue
            pool = stats.get(worker_name, {}).get("pool", {})
            concurrency = pool.get("max-concurrency")
            if concurrency:
                slots += concurrency
        if slots > 0:
            logging.info("annotation_worker_slots: live pool = %d", slots)
            return slots
        logging.info("annotation_worker_slots: no annotation_workers seen, using fallback %d", fallback)
    except Exception:
        log_traceback()
    return fallback
