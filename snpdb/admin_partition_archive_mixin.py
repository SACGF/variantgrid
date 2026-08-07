from django.contrib import messages

from snpdb.admin_utils import admin_action
from snpdb.partition_archive import (
    PartitionArchivePreconditionError,
    archive_partitioned_model,
)


class ArchivePartitionDataAdminMixin:
    """ Adds an "Archive partition data" admin action to any ModelAdmin
        whose model uses DataArchiveMixin + RelatedModelsPartitionModel.

        @see claude/issue_1537_archive_plan.md
    """

    @admin_action("Archive partition data")
    def archive_partition_data(self, request, queryset):
        for obj in queryset:
            try:
                archive = archive_partitioned_model(obj, request.user, reason="Admin action")
                messages.info(
                    request,
                    f"Queued archive for {obj} (PartitionArchive pk={archive.pk}). "
                    f"Ensure PARTITION_ARCHIVE_DIR has free space for the dump.",
                )
            except PartitionArchivePreconditionError as e:
                messages.error(request, f"{obj}: {e}")
