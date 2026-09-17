import logging
import operator
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta
from functools import reduce
from typing import Any

from django.conf import settings
from django.contrib import messages
from django.db import connection
from django.shortcuts import render
from django.utils.timesince import timesince
from django.utils.timezone import localtime

from annotation.models import (
    AnnotationRun,
    AnnotationStatus,
    VariantAnnotation,
    VariantAnnotationVersion,
)
from classification.models.classification_import_run import ClassificationImportRun
from genes.hgvs import HGVSMatcher
from library.django_utils import highest_pk, require_superuser
from library.health_check import HealthCheckRequest, health_check_overall_stats_signal
from library.integration_status import get_integration_statuses, run_integration_trigger
from library.log_utils import AdminNotificationBuilder, report_message, slack_bot_username
from library.utils import flatten_nested_lists
from snpdb.models import VCF, Sample, Variant
from snpdb.models.models_genome import GenomeBuild
from upload.upload_stats import get_vcf_variant_upload_stats
from variantgrid.celery import app
from variantgrid.tasks.server_monitoring_tasks import get_disk_messages
from variantopedia.server_status import get_dashboard_notices
from variantopedia.tasks.server_status_tasks import notify_server_status_now


def strip_celery_from_keys(celery_state):
    worker_status = {}
    if celery_state:
        for worker_string, data in celery_state.items():
            m = re.match(r".*@(.*?)$", worker_string)
            if m:
                worker = m.group(1)
                worker_status[worker] = data
    return worker_status


def _celery_worker_status() -> dict:
    if not settings.CELERY_ENABLED:
        return {}

    celery_workers = {}
    # This relies on the services being started with the "-n worker_name" with a separate one for each service
    worker_names = settings.CELERY_WORKER_NAMES.copy()
    if settings.URLS_APP_REGISTER["analysis"]:
        worker_names.extend(settings.CELERY_ANALYSIS_WORKER_NAMES)

    i = app.control.inspect()
    ping = strip_celery_from_keys(i.ping())
    stats = strip_celery_from_keys(i.stats())
    active = strip_celery_from_keys(i.active())
    scheduled = strip_celery_from_keys(i.scheduled())

    for worker in worker_names:
        num_workers = "?"
        status = 'ERROR - no workers found'
        ok = False

        # Sometimes stats fails - just use ping
        pong = ping.get(worker, {})
        if pong.get("ok") == "pong":
            status = "OK"
            ok = True

        if data := stats.get(worker):
            processes = data.get("pool", {}).get("processes")
            if processes:
                num_workers = len(processes)
                status = "OK"
                ok = True

        num_active = 0
        active_jobs = []
        if worker_active := active.get(worker):
            num_active = len(worker_active)
            logging.debug("worker %s active: %s", worker, worker_active)
            name_time_stamps = defaultdict(list)
            for worker_data in worker_active:
                name = worker_data.get("name")
                time_start = worker_data.get("time_start")
                if name and time_start:
                    name = name.split(".")[-1]  # remove prefix
                    name_time_stamps[name].append(time_start)

            for name, time_stamps in sorted(name_time_stamps.items(), key=lambda x: len(x[1]), reverse=True):
                earliest_ts = min(time_stamps)
                dt = datetime.fromtimestamp(earliest_ts)
                earliest = f"{timesince(dt)} ago"
                active_jobs.append(f"{name}: {len(time_stamps)} (earliest={earliest})")

        celery_workers[worker] = {
            "status": status,
            "ok": ok,
            "active": num_active,
            "scheduled": len(scheduled.get(worker, [])),
            "num_workers": num_workers,
            "active_jobs": ", ".join(active_jobs),
        }

    return celery_workers


@require_superuser
def server_status(request):
    if request.method == "POST":
        action = request.POST.get('action')
        if action == 'Test Slack':
            nb = AdminNotificationBuilder(message="Slack Check")
            nb.add_markdown("This is a Slack Test :ladybug:")
            nb.send()
            messages.add_message(request, level=messages.INFO, message="Slack should have been sent a test message.")
        elif action == 'Health Check':
            notify_server_status_now()
            messages.add_message(request, level=messages.INFO, message="Slack should have been sent the health check.")
        elif action == 'Test Rollbar':
            report_message("Testing Rollbar", level='error')
            messages.add_message(request, level=messages.INFO, message="Rollbar should have been sent an error.")
        elif action == 'Test Message Branding':
            messages.success(request, "Success message")
            messages.info(request, "Info message")
            messages.warning(request, "Warning message")
            messages.error(request, "Error message")

        elif action == 'run-integration':
            action_id = request.POST.get('action_id')
            if matched := run_integration_trigger(action_id):
                integration, trigger = matched
                messages.add_message(request, level=messages.INFO,
                                     message=f"Ran '{trigger.label}' for {integration.name}.")
            else:
                messages.add_message(request, level=messages.ERROR,
                                     message=f"No integration is registered against '{action_id}'.")
        elif action == 'kill-pid':
            pid = int(request.POST.get('pid'))
            with connection.cursor() as cursor:
                cursor.execute("SELECT pg_terminate_backend(%s)", [pid])
                terminated = cursor.fetchone()[0]
                messages.add_message(request, level=messages.INFO, message=f"Query {pid} Terminated = {terminated}")
        else:
            logging.warning("Unrecognised action %s", action)

        # return redirect(reverse('server_status'))

        # TODO should redirect to read-only version of the page

    celery_workers = _celery_worker_status()

    can_access_reference = True
    try:
        for genome_build in GenomeBuild.builds_with_annotation():
            _ = genome_build.reference_fasta  # Throws exception on error
    except (KeyError, FileNotFoundError):
        can_access_reference = False

    # Variant Annotation - incredibly quick check
    highest_variant_annotated = {}
    try:
        q = reduce(operator.and_, VariantAnnotation.VARIANT_ANNOTATION_Q)
        highest_variant = Variant.objects.filter(q).order_by("pk").last()
        genome_build = next(iter(highest_variant.genome_builds))  # Just pick one if spans multiple
        vav = VariantAnnotationVersion.latest(genome_build)
        annotated = highest_variant.variantannotation_set.filter(version=vav).exists()
        if annotated:
            highest_variant_annotated["status"] = "info"
            highest_variant_annotated["message"] = "OK"
        else:
            try:
                ar_qs = AnnotationRun.objects.filter(annotation_range_lock__version=vav,
                                                     annotation_range_lock__min_variant__lte=highest_variant.pk,
                                                     annotation_range_lock__max_variant__gte=highest_variant.pk)
                ar_qs = ar_qs.exclude(status=AnnotationStatus.ERROR)
                annotation_run_message = f"AnnotationRuns: {', '.join([str(ar) for ar in ar_qs])}"

                highest_variant_annotated["status"] = "warning"
                highest_variant_annotated["message"] = annotation_run_message
            except AnnotationRun.DoesNotExist:
                highest_variant_annotated["status"] = "danger"
                highest_variant_annotated["message"] = "Not annotated, no AnnotationRun!"
    except Exception as e:
        highest_variant_annotated["status"] = "danger"
        highest_variant_annotated["message"] = str(e)

    disk_messages = get_disk_messages(info_messages=True)
    disk_free = {"status": "info", "messages": []}
    for status, message in disk_messages:
        if status == "warning":
            disk_free["status"] = "warning"
        disk_free["messages"].append(message)

    integrations = []
    integration_messages = []
    if settings.INTEGRATION_STATUS_ENABLED:
        integrations, integration_messages = get_integration_statuses()

    context = {
        "celery_workers": celery_workers,
        "queries": long_running_sql(0),
        "can_access_reference": can_access_reference,
        "highest_variant_annotated": highest_variant_annotated,
        "integrations": integrations,
        "integration_messages": integration_messages,
    }
    return render(request, "variantopedia/server_status.html", context)


@require_superuser
def server_status_activity(request, days_ago: int):
    dashboard_notices = get_dashboard_notices(request.user, days_ago)
    return render(request, "variantopedia/server_status_activity_detail.html", {"dashboard_notices": dashboard_notices})


@require_superuser
def server_status_settings(request):
    slack_emoji = (settings.SLACK or {}).get('emoji') or ':dna:'
    slack_username = f"{slack_emoji} {slack_bot_username()}"

    hgvs_matcher = HGVSMatcher(GenomeBuild.grch38())

    return render(request, "variantopedia/server_status_settings_detail.html", {
        "settings": settings,
        "slack_bot_username": slack_username,
        "ongoing_imports": ClassificationImportRun.ongoing_imports(),
        "hgvs_matcher": hgvs_matcher
    })


@require_superuser
def health_check_details(request):
    now = localtime()
    since = now - timedelta(days=1)
    health_request = HealthCheckRequest(since=since, now=now)

    results = []
    for _, result in health_check_overall_stats_signal.send_robust(sender=None, health_request=health_request):
        if not isinstance(result, Exception):
            results.append(result)

    checks = flatten_nested_lists(results)
    checks = sorted(checks, key=lambda hc: (hc.sort_order(), hc.name if hasattr(hc, "name") else "z"))
    overall_lines = []
    for check in checks:
        line_content = check.as_html()
        overall_lines.append(line_content)

    context = {
        'overall_lines': overall_lines,
    }
    return render(request, "variantopedia/health_check_details.html", context)


@dataclass
class RunningQuery:
    pid: int
    duration: Any  # is actually a Duration
    query: str
    state: str


def long_running_sql(min_age_in_seconds: int = 30):
    with connection.cursor() as cursor:
        db_name = connection.settings_dict['NAME']
        # We exclude "idle" state - as they are from connection pool and not actually running
        cursor.execute(
            """
            SELECT
              pid,
              now() - pg_stat_activity.query_start AS duration,
              query,
              state
            FROM pg_stat_activity
            WHERE (now() - pg_stat_activity.query_start) > interval %s
            AND datname = %s
            AND state <> 'idle'                
            ORDER BY now() - pg_stat_activity.query_start desc;
            """,
            [f"{min_age_in_seconds} seconds", db_name]
        )

        def to_obj(row) -> RunningQuery:
            return RunningQuery(
                pid=row[0],
                duration=row[1],
                query=row[2],
                state=row[3]
            )

        return [to_obj(result) for result in cursor.fetchall()]


@require_superuser
def database_statistics(request):
    max_variant_id = highest_pk(Variant)
    num_vcfs = VCF.objects.count()
    num_samples = Sample.objects.count()

    variant_stats_per_build = defaultdict(dict)
    for genome_build in GenomeBuild.builds_with_annotation():
        vcf_variant_stats_df = get_vcf_variant_upload_stats(genome_build)
        for col in ["cumulative_samples", "cumulative_variants", "cumulative_genotypes", "percent_known"]:
            variant_stats_per_build[genome_build.name][col] = vcf_variant_stats_df[col].tolist()

    context = {"max_variant_id": max_variant_id,
               "num_vcfs": num_vcfs,
               "num_samples": num_samples,
               "variant_stats_per_build": dict(variant_stats_per_build)}
    return render(request, "variantopedia/database_statistics_detail.html", context)
