"""
Whether this box survives running out of memory: Redis protected from the kernel OOM killer (it is often the
largest single process, and losing it fails every in-flight task at once), and swap to absorb a spike (#2187).

Entry point: `check_memory_protection`, run by `manage.py deployment_check`.
"""
import os

import redis
from django.conf import settings

REDIS_OOM_DROP_IN = "config/systemd/redis-server.service.d/variantgrid.conf"


def _redis_oom_score_adj() -> int | None:
    """ None when Redis isn't on this box (its pid isn't ours to read) """
    client = redis.Redis.from_url(settings.CACHES["default"]["LOCATION"])
    redis_pid = client.info("server")["process_id"]
    try:
        with open(f"/proc/{redis_pid}/oom_score_adj") as f:
            return int(f.read())
    except FileNotFoundError:
        return None


def _swap_total_kb() -> int:
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("SwapTotal:"):
                return int(line.split()[1])
    return 0


def check_memory_protection() -> dict:
    if not os.path.exists("/proc/meminfo"):
        return {}  # not Linux

    # Warnings, not invalid: the upgrader dies on an invalid check, and neither stops the box working
    redis_check = {"valid": True}
    oom_score_adj = _redis_oom_score_adj()
    if oom_score_adj is not None and oom_score_adj >= 0:
        redis_check["warning"] = (f"Redis oom_score_adj={oom_score_adj} - the OOM killer may pick it before a worker, "
                                  f"install {REDIS_OOM_DROP_IN}")

    swap_check = {"valid": True}
    if not _swap_total_kb():
        swap_check["warning"] = "No swap - a memory spike (eg a sequencing run's imports) goes straight to the OOM killer"

    return {
        "Redis OOM score": redis_check,
        "Swap": swap_check,
    }
