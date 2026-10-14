"""
The systemd units a deployment runs: gunicorn, celery beat and one celeryd_<queue> worker unit per celery
queue (`config/systemd/`). Shared by `vg status` and the deployment checks that read the installed units.
"""
from django.conf import settings

SERVICES = ("gunicorn", "celeryd_beat")


def service_names() -> list[str]:
    return [*SERVICES, *(f"celeryd_{queue.name}" for queue in settings.CELERY_TASK_QUEUES)]
