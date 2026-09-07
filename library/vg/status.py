"""
`vg status`: the running deployment in one call - settings module, git, database size and the big
tables, current annotation version per build, celery services and queue depths, annotation runs in
flight, recent ERROR events, outstanding manual migration tasks and free disk on the data roots.

Each section is gathered independently and a failure becomes that section's `error` rather than
aborting the whole report, because the point of the command is to work when the box is unwell.
Read-only: nothing here writes to the database, the broker or the filesystem.
"""
import shutil
import socket
import subprocess
from dataclasses import dataclass, field
from typing import Any

from django.conf import settings
from django.db import connection
from django.db.models import Count
from django.utils import timezone

from annotation.models import AnnotationRun, VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus
from eventlog.models import Event
from library.enums.log_level import LogLevel
from library.vg.repo import REPO_ROOT, git
from library.vg.settings_chain import resolved_settings_module
from manual.models import ManualMigrationOutstanding
from snpdb.models import GenomeBuild
from variantgrid.celery import app as celery_app

# pg_class estimates for the tables whose size shapes how a change must be written (see claude/guides/operations.md#scale)
SCALE_TABLES = ("snpdb_variant", "snpdb_allele", "snpdb_locus", "snpdb_sample", "snpdb_vcf",
                "classification_classification", "analysis_analysis", "analysis_varianttag", "annotation_annotationrun")
SERVICES = ("gunicorn", "celeryd_beat")
RECENT_ERRORS = 10


@dataclass
class Status:
    host: str
    settings_module: str
    git: dict[str, Any]
    db: dict[str, Any] = field(default_factory=dict)
    builds: dict[str, Any] = field(default_factory=dict)
    services: dict[str, str] = field(default_factory=dict)
    queues: dict[str, Any] = field(default_factory=dict)
    annotation_runs: dict[str, int] = field(default_factory=dict)
    manual_outstanding: list[str] = field(default_factory=list)
    recent_errors: list[dict[str, str]] = field(default_factory=list)
    disk: dict[str, str] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)


def _section(status: Status, name: str, gather):
    try:
        setattr(status, name, gather())
    except Exception as e:  # pylint: disable=broad-exception-caught
        status.errors[name] = f"{type(e).__name__}: {e}"


def gather_status() -> Status:
    status = Status(host=socket.gethostname(), settings_module=resolved_settings_module(), git=_git())
    _section(status, "db", _db)
    _section(status, "builds", _builds)
    _section(status, "services", _services)
    _section(status, "queues", _queues)
    _section(status, "annotation_runs", _annotation_runs)
    _section(status, "manual_outstanding", _manual_outstanding)
    _section(status, "recent_errors", _recent_errors)
    _section(status, "disk", _disk)
    return status


def _git() -> dict[str, Any]:
    try:
        return {"sha": git("rev-parse", "--short", "HEAD").strip(),
                "branch": git("rev-parse", "--abbrev-ref", "HEAD").strip(),
                "dirty": bool(git("status", "--porcelain").strip())}
    except (OSError, subprocess.CalledProcessError) as e:
        return {"error": str(e)}


def _db() -> dict[str, Any]:
    with connection.cursor() as cursor:
        cursor.execute("SELECT current_database(), pg_database_size(current_database()), version()")
        name, size_bytes, version = cursor.fetchone()
        cursor.execute("SELECT relname, reltuples::bigint FROM pg_class WHERE relname = ANY(%s)", [list(SCALE_TABLES)])
        estimates = dict(cursor.fetchall())
    counts = {table: estimates[table] for table in SCALE_TABLES if table in estimates}
    return {"name": name, "size_gb": round(size_bytes / 1e9, 1), "postgres": version.split(",")[0],
            "estimated_rows": counts}


def _builds() -> dict[str, Any]:
    builds = {}
    for genome_build in GenomeBuild.builds_with_annotation().order_by("name"):
        vav = VariantAnnotationVersion.latest(genome_build)
        builds[genome_build.name] = {
            "vav": vav.pk, "vep": vav.vep, "columns_version": vav.columns_version,
            "annotation_date": vav.annotation_date.date().isoformat(),
        } if vav else {"vav": None}
    return builds


def _service_names() -> list[str]:
    return [*SERVICES, *(f"celeryd_{queue.name}" for queue in settings.CELERY_TASK_QUEUES)]


def _services() -> dict[str, str]:
    if shutil.which("systemctl") is None:
        return {"error": "systemctl not available"}
    names = _service_names()
    result = subprocess.run(["systemctl", "is-active", *names], capture_output=True, text=True, check=False)
    states = result.stdout.split()
    return dict(zip(names, states)) if len(states) == len(names) else {"error": result.stderr.strip() or result.stdout.strip()}


def _queues() -> dict[str, Any]:
    """ Message and consumer counts per queue via a passive declare on the broker (no queue is created) """
    queues = {}
    with celery_app.connection_for_read() as connection_, connection_.channel() as channel:
        for queue in settings.CELERY_TASK_QUEUES:
            try:
                declared = channel.queue_declare(queue=queue.name, passive=True)
                queues[queue.name] = {"messages": declared.message_count, "consumers": declared.consumer_count}
            except Exception as e:  # pylint: disable=broad-exception-caught
                queues[queue.name] = {"error": type(e).__name__}
    return queues


def _annotation_runs() -> dict[str, int]:
    in_flight = AnnotationRun.objects.exclude(status__in=AnnotationStatus.get_completed_states())
    by_status = {}
    for run_status, count in in_flight.values_list("status").annotate(n=Count("pk")).order_by("status"):
        by_status[AnnotationStatus(run_status).label] = count
    return by_status


def _manual_outstanding() -> list[str]:
    return [task.to_json()["line"] for task in ManualMigrationOutstanding.outstanding_tasks()]


def _recent_errors() -> list[dict[str, str]]:
    events = Event.objects.filter(severity=LogLevel.ERROR).order_by("-date")[:RECENT_ERRORS]
    return [{"date": timezone.localtime(e.date).strftime("%Y-%m-%d %H:%M"), "app": e.app_name, "name": e.name,
             "details": (e.details or "").strip().replace("\n", " ")[:160]} for e in events]


def _disk() -> dict[str, str]:
    roots = {settings.PRIVATE_DATA_ROOT, settings.ANNOTATION_BASE_DIR, settings.MEDIA_ROOT, str(REPO_ROOT)}
    disk = {}
    for root in sorted(roots):
        try:
            usage = shutil.disk_usage(root)
            disk[root] = f"{usage.free / 1e9:.0f} GB free of {usage.total / 1e9:.0f} GB ({usage.used / usage.total:.0%} used)"
        except OSError as e:
            disk[root] = f"error: {e.strerror}"
    return disk


def render_status(status: Status) -> str:
    lines = [f"{status.host}  settings={status.settings_module}  git={status.git.get('branch')}@{status.git.get('sha')}"
             + ("  (dirty)" if status.git.get("dirty") else "")]
    if status.db:
        db = status.db
        lines.append(f"db: {db.get('name')} {db.get('size_gb')} GB  {db.get('postgres')}")
        rows = db.get("estimated_rows", {})
        if rows:
            lines.append("  rows≈ " + ", ".join(f"{table.split('_', 1)[1]} {n:,}" for table, n in rows.items()))
    for build, info in status.builds.items():
        lines.append(f"{build}: VAV {info.get('vav')}  VEP {info.get('vep')}  columns v{info.get('columns_version')}  "
                     f"{info.get('annotation_date')}")
    if status.services:
        down = [name for name, state in status.services.items() if state != "active"]
        lines.append(f"services: {len(status.services) - len(down)}/{len(status.services)} active"
                     + (f"  DOWN: {', '.join(down)}" if down else ""))
    if status.queues:
        lines.append("queues: " + ", ".join(
            f"{name} {info['messages']}q/{info['consumers']}c" if "messages" in info else f"{name} {info['error']}"
            for name, info in status.queues.items()))
    lines.append("annotation runs in flight: " + (", ".join(f"{k} {v}" for k, v in status.annotation_runs.items()) or "none"))
    lines.append(f"manual tasks outstanding: {len(status.manual_outstanding)}"
                 + ("".join(f"\n  {line}" for line in status.manual_outstanding) if status.manual_outstanding else ""))
    lines.append(f"recent ERROR events ({len(status.recent_errors)}):")
    lines += [f"  {e['date']} {e['app']}/{e['name']}: {e['details']}" for e in status.recent_errors]
    for root, usage in status.disk.items():
        lines.append(f"disk {root}: {usage}")
    for name, error in status.errors.items():
        lines.append(f"!! {name}: {error}")
    return "\n".join(lines)
