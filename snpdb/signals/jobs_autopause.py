import logging
from typing import Optional

from celery.signals import worker_ready
from django.conf import settings

from library.log_utils import log_traceback
from snpdb.models import JobsControl


def _read_proc_uptime_seconds() -> float:
    with open("/proc/uptime") as f:
        return float(f.read().split()[0])


def _read_proc_boot_time() -> Optional[int]:
    with open("/proc/stat") as f:
        for line in f:
            if line.startswith("btime "):
                return int(line.split()[1])
    return None


@worker_ready.connect
def autopause_jobs_after_reboot(**kwargs):
    """ Crash safety brake. On a long-lived host a low system uptime means the box rebooted - an
        ordinary deploy restarts the worker process but leaves uptime high. After a reboot we pause
        the analysis + annotation dispatchers once per boot (keyed on /proc/stat btime so concurrent
        workers and later restarts don't re-trip it), so jobs that may have crashed the machine don't
        immediately re-launch. An admin resumes via 'manage.py jobs_control resume'.

        Disable on ephemeral / autoscaled hosts (where a fresh boot is routine) via
        JOBS_AUTOPAUSE_ON_REBOOT. """
    if not getattr(settings, "JOBS_AUTOPAUSE_ON_REBOOT", False):
        return
    try:
        uptime = _read_proc_uptime_seconds()
        if uptime >= settings.JOBS_AUTOPAUSE_ON_REBOOT_UPTIME_SECS:
            return  # box up a while - ordinary worker / deploy restart, not a reboot
        boot_time = _read_proc_boot_time()
        if boot_time is None:
            return
        reason = (f"Auto-paused: host booted {int(uptime)}s ago (possible crash recovery). "
                  f"Inspect, then run 'manage.py jobs_control resume'.")
        if JobsControl.autopause_for_boot(boot_time, reason=reason, by="auto:reboot"):
            logging.warning("JobsControl %s", reason)
    except Exception:
        log_traceback()  # never block worker startup
