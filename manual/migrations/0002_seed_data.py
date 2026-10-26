from django.db import migrations

from manual.operations.manual_operations import ManualOperation
from variantgrid.deployment_validation.memory_check import gunicorn_unit_missing_memory_cap, redis_unprotected_from_oom_killer
from variantgrid.deployment_validation.systemd_check import systemd_env_files_problem


def _gunicorn_unit_needs_memory_cap(_apps) -> bool:
    """ Only a box whose installed gunicorn unit lacks the cap has anything to do: a dev box with no unit,
        or a deployment that has already copied the new file, registers nothing """
    return gunicorn_unit_missing_memory_cap() is not None


def _redis_needs_oom_drop_in(_apps) -> bool:
    """ Only a box running Redis itself with the default OOM score has anything to do: a deployment with a
        remote Redis, or one that has already installed the drop-in, registers nothing """
    return redis_unprotected_from_oom_killer()


def _systemd_env_files_needed(_apps) -> bool:
    """ Only a box whose installed units read an env file that isn't there, or whose VG_INSTALL_DIR is another
        checkout, has anything to do: a dev box with no units, or a deployment already set up, registers nothing """
    return systemd_env_files_problem() is not None


class Migration(migrations.Migration):

    initial = True  # manual.squashed_history records it with the initials on a reset database

    dependencies = [
        ('manual', '0001_initial'),
    ]

    operations = [
        ManualOperation.operation_other(
            args=["Install config/systemd/gunicorn.service (MemoryHigh/MemoryMax/OOMPolicy=continue): "
                  "sudo cp config/systemd/gunicorn.service /lib/systemd/system/ && "
                  "sudo systemctl daemon-reload && sudo systemctl restart gunicorn. The unit reads "
                  "/etc/variantgrid/variantgrid.env - if that is missing, do the env files step first "
                  "(config/systemd/README.txt)"],
            note="Caps the gunicorn unit's memory at a percentage of the box's RAM so a request that pulls a whole "
                 "table into a worker loses that worker, not the box. deployment_check warns until it is installed.",
            test=_gunicorn_unit_needs_memory_cap),
        ManualOperation.operation_other(
            args=["Install config/systemd/redis-server.service.d/variantgrid.conf (OOMScoreAdjust=-900) as a "
                  "drop-in for the distribution's redis unit - commands per distribution in "
                  "config/systemd/README.txt - then daemon-reload and restart redis"],
            note="Redis is often the largest single process, so the kernel OOM killer picks it first and every "
                 "in-flight task loses its cache, locks and results at once (#2187). The drop-in makes it pick a "
                 "celery worker instead. deployment_check warns until it is installed.",
            test=_redis_needs_oom_drop_in),
        ManualOperation.operation_other(
            args=["Install the systemd env files the units read: "
                  "sudo cp -rn config/celery config/gunicorn config/variantgrid.env /etc/variantgrid/ then set "
                  "VG_INSTALL_DIR in /etc/variantgrid/variantgrid.env to this checkout, "
                  "sudo cp config/systemd/*.service /lib/systemd/system/ && sudo systemctl daemon-reload, "
                  "then restart gunicorn and the celeryd_* units"],
            note="The units in config/systemd are one file for every install location: they get the install dir from "
                 "/etc/variantgrid/variantgrid.env and the gunicorn/celery binaries from their own .env. Installed "
                 "without those files a unit fails to start. deployment_check is invalid until they are in place.",
            test=_systemd_env_files_needed),
    ]
