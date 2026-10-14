"""
Whether this box survives running out of memory: Redis protected from the kernel OOM killer (it is often the
largest single process, and losing it fails every in-flight task at once), swap to absorb a spike (#2187), and
the gunicorn unit capped so one request that pulls a whole table into a worker loses that worker, not the box.

Entry point: `check_memory_protection`, run by `manage.py deployment_check`. `redis_unprotected_from_oom_killer`
and `gunicorn_unit_missing_memory_cap` are shared with the manual migrations that register the deploy steps.
"""
import os

import redis
from django.conf import settings

from variantgrid.deployment_validation.systemd_check import systemctl_show

REDIS_OOM_DROP_IN = "config/systemd/redis-server.service.d/variantgrid.conf"
GUNICORN_UNIT = "config/systemd/gunicorn.service"


def _redis_oom_score_adj() -> int | None:
    """ None when Redis isn't on this box (its pid isn't ours to read) """
    client = redis.Redis.from_url(settings.CACHES["default"]["LOCATION"])
    redis_pid = client.info("server")["process_id"]
    try:
        with open(f"/proc/{redis_pid}/oom_score_adj", encoding="utf-8") as f:
            return int(f.read())
    except FileNotFoundError:
        return None


def redis_unprotected_from_oom_killer() -> bool:
    """ Redis runs on this box with a non-negative OOM score, so the kernel may pick it before a worker.
        False when Redis is remote or unreachable, or the cache isn't Redis at all (the test settings):
        none is a reason to fail a migrate """
    try:
        oom_score_adj = _redis_oom_score_adj()
    except (redis.RedisError, ValueError):
        return False
    return oom_score_adj is not None and oom_score_adj >= 0


def _swap_total_kb() -> int:
    with open("/proc/meminfo", encoding="utf-8") as f:
        for line in f:
            if line.startswith("SwapTotal:"):
                return int(line.split()[1])
    return 0


def gunicorn_unit_memory_settings() -> dict[str, str] | None:
    """ MemoryMax / OOMPolicy of the installed gunicorn systemd unit, or None where there is no such unit
        (a dev box running runserver, a host without systemd) and so nothing to check """
    props = systemctl_show("gunicorn", "MemoryMax", "OOMPolicy")
    if props is None:
        return None
    return {key: values[0] for key, values in props.items()}


def gunicorn_unit_missing_memory_cap() -> str | None:
    """ What the installed gunicorn unit lacks, or None when it is capped as the repo unit is """
    props = gunicorn_unit_memory_settings()
    if props is None:
        return None
    missing = []
    if props.get("MemoryMax") == "infinity":
        missing.append("no MemoryMax")
    if props.get("OOMPolicy") != "continue":
        missing.append(f"OOMPolicy={props.get('OOMPolicy')}")
    return ", ".join(missing) or None


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

    gunicorn_check = {"valid": True}
    if missing := gunicorn_unit_missing_memory_cap():
        gunicorn_check["warning"] = (f"gunicorn unit has {missing} - one request loading a big table takes the box "
                                     f"down rather than its worker, install {GUNICORN_UNIT}")

    return {
        "Redis OOM score": redis_check,
        "Swap": swap_check,
        "gunicorn memory cap": gunicorn_check,
    }
