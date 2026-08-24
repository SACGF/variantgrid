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
    """ Crash safety brake. The worker records each host boot it starts under (/proc/stat btime),
        and pauses the analysis + annotation dispatchers when the previous boot was recent - a box
        that comes up and goes down again inside the window is likely being crashed by the jobs it
        relaunches. A single reboot the box survives leaves the dispatchers running, so an ordinary
        reboot auto-starts. A paused deployment is resumed via 'manage.py jobs_control resume'.

        On a long-lived host a low system uptime means the box rebooted - an ordinary deploy
        restarts the worker process but leaves uptime high - so a worker that starts on a box which
        has been up a while records the boot and leaves the dispatchers alone.

        Disable on ephemeral / autoscaled hosts (where a fresh boot is routine) via
        JOBS_AUTOPAUSE_ON_REBOOT. """
    if not getattr(settings, "JOBS_AUTOPAUSE_ON_REBOOT", False):
        return
    try:
        boot_time = _read_proc_boot_time()
        if boot_time is None:
            return
        registration = JobsControl.register_boot(boot_time)
        if not registration.new_boot:
            return  # already judged this boot - a concurrent worker or a later worker restart
        if registration.previous_boot_time is None:
            return  # first boot we've seen - nothing to measure a crash loop against
        uptime = _read_proc_uptime_seconds()
        if uptime >= settings.JOBS_AUTOPAUSE_ON_REBOOT_UPTIME_SECS:
            return  # box up a while - ordinary worker / deploy restart, not a reboot
        secs_since_previous_boot = boot_time - registration.previous_boot_time
        if secs_since_previous_boot > settings.JOBS_AUTOPAUSE_ON_REBOOT_WINDOW_SECS:
            return  # last boot was long ago - a one-off reboot, let the jobs start
        reason = (f"Auto-paused: host rebooted twice within "
                  f"{secs_since_previous_boot}s (possible crash loop). "
                  f"Inspect, then run 'manage.py jobs_control resume'.")
        JobsControl.pause(reason=reason, by="auto:reboot")
        logging.warning("JobsControl %s", reason)
    except Exception:
        log_traceback()  # never block worker startup
