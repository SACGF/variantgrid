"""
Celery task that performs the long-running pg_dump+drop for a PartitionArchive row.

@see claude/issue_1537_archive_plan.md
"""

import hashlib
import logging
import os
import shlex
import subprocess

import celery
from django.conf import settings
from django.db import transaction
from django.utils import timezone

from eventlog.models import create_event
from library.utils.database_utils import run_sql
from snpdb.models.models_partition_archive import PartitionArchive


@celery.shared_task(queue="db_workers")
def perform_partition_archive(archive_pk: int, reason: str = ""):
    archive = PartitionArchive.objects.get(pk=archive_pk)
    if archive.status != PartitionArchive.Status.PENDING:
        logging.warning(
            "perform_partition_archive: archive %s is %s, not PENDING -- skipping",
            archive_pk, archive.status,
        )
        return

    archive.status = PartitionArchive.Status.IN_PROGRESS
    archive.started_at = timezone.now()
    archive.save(update_fields=["status", "started_at"])

    source = archive.resolve_source()
    if source is None:
        _fail(archive, f"Source {archive.source_label} not found")
        return

    child_tables = [source.get_partition_table(base_table_name=t)
                    for t in archive.source_table_names]

    try:
        _run_pg_dump(archive.dump_path, child_tables)
        _verify_dump(archive.dump_path, child_tables)
        archive.sha256 = _sha256_file(archive.dump_path)
        archive.save(update_fields=["sha256"])

        with transaction.atomic():
            for child in child_tables:
                run_sql(f'DROP TABLE IF EXISTS "{child}";')

            source.data_archived_date = timezone.now()
            source.data_archived_by = archive.dumped_by
            source.data_archive_reason = reason
            source.data_restorable_from = archive.dump_path
            source.save(update_fields=["data_archived_date", "data_archived_by",
                                       "data_archive_reason", "data_restorable_from"])

            archive.status = PartitionArchive.Status.COMPLETE
            archive.completed_at = timezone.now()
            archive.save(update_fields=["status", "completed_at"])

        create_event(
            archive.dumped_by, "partition_archive_complete",
            details=(f"archive_pk={archive.pk} archive_name={archive.archive_name!r} "
                     f"source={archive.source_label} dump_path={archive.dump_path!r} "
                     f"sha256={archive.sha256} reason={reason!r}"),
            app_name="snpdb",
        )
    except Exception as exc:
        logging.exception("perform_partition_archive failed for %s", archive_pk)
        _fail(archive, str(exc))


def _run_pg_dump(dump_path: str, tables: list[str]):
    db = settings.DATABASES["default"]
    table_args = []
    for t in tables:
        table_args += ["--table", t]
    cmd = [
        "pg_dump",
        "--host", db.get("HOST") or "localhost",
        "--port", str(db.get("PORT") or 5432),
        "--username", db["USER"],
        "--dbname", db["NAME"],
        "--format", "custom",
        "--file", dump_path,
        *table_args,
    ]
    env = os.environ.copy()
    if db.get("PASSWORD"):
        env["PGPASSWORD"] = db["PASSWORD"]
    logging.info("pg_dump: %s", " ".join(shlex.quote(c) for c in cmd))
    subprocess.run(cmd, check=True, env=env, capture_output=True, text=True)


def _verify_dump(dump_path: str, expected_tables: list[str]):
    """ Run `pg_restore --list` and confirm every expected table appears. """
    result = subprocess.run(["pg_restore", "--list", dump_path],
                            check=True, capture_output=True, text=True)
    listing = result.stdout
    missing = [t for t in expected_tables if t not in listing]
    if missing:
        raise RuntimeError(f"pg_restore --list missing tables: {missing}")


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _fail(archive: PartitionArchive, msg: str):
    archive.status = PartitionArchive.Status.FAILED
    archive.error_message = msg
    archive.completed_at = timezone.now()
    archive.save(update_fields=["status", "error_message", "completed_at"])
    create_event(archive.dumped_by, "partition_archive_failed",
                 details=f"archive_pk={archive.pk} error={msg!r}",
                 app_name="snpdb")
