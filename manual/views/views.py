from typing import Dict, Any

from django.conf import settings
from django.contrib import messages
from django.contrib.messages import add_message
from django.db.models import QuerySet
from django.http import HttpRequest, HttpResponseRedirect
from django.shortcuts import render
from django.urls import reverse

from library.django_utils import require_superuser
from library.git import Git
from manual.models import ManualMigrationAttempt, ManualMigrationTask, ManualMigrationOutstanding
from snpdb.views.datatable_view import DatatableConfig, RichColumn, SortOrder

"""
    task = models.ForeignKey(ManualMigrationTask, on_delete=CASCADE)
    source_version = models.TextField(null=True, blank=True)
    requires_retry = models.BooleanField(default=False)
    note = models.TextField(null=True, blank=True)
"""


class MigrationAttemptColumns(DatatableConfig[ManualMigrationAttempt]):

    def render_task(self, row: Dict[str, Any]) -> str:
        task_id = row["task_id"]
        return ManualMigrationTask.describe_manual(task_id)

    def __init__(self, request: HttpRequest):
        super().__init__(request)

        self.rich_columns = [
            RichColumn(
                key="task_id",
                name="task",
                orderable=True,
                renderer=self.render_task,
                css_class="text-monospace"
            ),
            RichColumn(key="source_version", label="Git Version", orderable=True),
            RichColumn(key="note", orderable=True),
            RichColumn(key="created", client_renderer='TableFormat.timestamp', orderable=True, default_sort=SortOrder.DESC)
        ]

    def get_initial_queryset(self) -> QuerySet[ManualMigrationAttempt]:
        return ManualMigrationAttempt.objects.all()

    def filter_queryset(self, qs: QuerySet[ManualMigrationAttempt]) -> QuerySet[ManualMigrationAttempt]:
        if self.get_query_param("exclude_standard") == "true":
            qs = qs.exclude(task_id__in=["git*pull", "manage*migrate", "manage*collectstatic"])
        return qs


@require_superuser
def manual_migrations_view(request: HttpRequest):
    outstanding_tasks = ManualMigrationOutstanding.outstanding_tasks()
    if request.method == "POST":
        if request.POST.get("action") == "skip":
            git = Git(settings.BASE_DIR)
            for out_task in outstanding_tasks:
                ManualMigrationAttempt.objects.create(
                    task=out_task.task,
                    source_version=git.hash,
                    requires_retry=False,
                    note=f"Manually skipped by {request.user.username}"
                )

            add_message(request, messages.SUCCESS, "Outstanding migration tasks marked as skipped")
            return HttpResponseRedirect(reverse('manual_migrations'))

    return render(request, "manual/manual_migrations.html", {
        "outstanding": outstanding_tasks
    })
