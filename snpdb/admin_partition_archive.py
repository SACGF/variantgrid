from django.contrib import admin, messages
from django.db.models import QuerySet

from snpdb.admin_utils import ModelAdminBasics, admin_action
from snpdb.models.models_partition_archive import PartitionArchive
from snpdb.partition_archive import (
    clear_partition_archive_dump,
    mark_partition_archive_restored,
)


@admin.register(PartitionArchive)
class PartitionArchiveAdmin(ModelAdminBasics):
    list_display = ("pk", "archive_name", "status", "source_app_label",
                    "source_model", "source_pk", "dumped_by", "dumped_date",
                    "completed_at", "restored_date", "cleared_date")
    list_filter = ("status", "source_app_label", "source_model")
    search_fields = ("archive_name", "source_pk")
    readonly_fields = ("archive_name", "dump_path", "sha256",
                       "source_app_label", "source_model", "source_pk",
                       "source_table_names", "status", "error_message",
                       "started_at", "completed_at", "dumped_by", "dumped_date",
                       "cleared_date", "cleared_by", "restored_date", "restored_by")

    def has_add_permission(self, request):
        return False

    @admin_action("Mark restored (after manual pg_restore)")
    def mark_restored(self, request, queryset: QuerySet[PartitionArchive]):
        for archive in queryset:
            try:
                mark_partition_archive_restored(archive, request.user)
                messages.info(request, f"Marked {archive} as restored")
            except ValueError as e:
                messages.error(request, f"{archive}: {e}")

    @admin_action("Clear dump file")
    def clear_dump(self, request, queryset: QuerySet[PartitionArchive]):
        for archive in queryset:
            try:
                clear_partition_archive_dump(archive, request.user)
                messages.info(request, f"Cleared dump for {archive}")
            except ValueError as e:
                messages.error(request, f"{archive}: {e}")
