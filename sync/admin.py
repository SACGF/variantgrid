import json
from typing import Optional

from django.contrib import admin, messages
from django.contrib.admin import RelatedFieldListFilter
from django.db.models import QuerySet

from snpdb.admin_utils import ModelAdminBasics, admin_action, admin_list_column
from sync.models import SyncRun, ClassificationModificationSyncRecord
from sync.models.models import SyncDestination
from sync.sync_runner import sync_runner_for_destination


@admin.register(SyncDestination)
class SyncDestinationAdmin(ModelAdminBasics):
    list_display = ('name', 'config', 'enabled')

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget()
        }, **kwargs)

    def _run_sync(self, request, queryset: QuerySet[SyncDestination], max_rows: Optional[int] = None):
        for sync_destination in queryset:
            sync_destination.run(full_sync=False, max_rows=max_rows)
            self.message_user(request, message=f"Completed {str(sync_destination)} row limit = {max_rows}")

    @admin_action("Validate configuration")
    def validate_configuration(self, request, queryset: QuerySet[SyncDestination]):
        for destination in queryset:
            if key := destination.config.get("sync_details"):
                try:
                    sync_details = destination.sync_details
                    sync_detail_keys = ", ".join(f"\"{key}\"" for key in sync_details.keys())
                    self.message_user(request, f"SyncDetails for \"{key}\" has entries for {sync_detail_keys} in the secrets file", level=messages.INFO)
                except ValueError:
                    self.message_user(request, f"SyncDetails for \"{key}\" are missing from the secrets file", level=messages.ERROR)
            else:
                self.message_user(request, f"SyncDestination \"{destination.name}\" has not specified a sync_details", level=messages.ERROR)

    # @admin_action("Row Counts")
    # def eligible_row_count(self, request, queryset):
    #     for destination in queryset:
    #         uploader = VariantGridUploadSyncer(destination)
    #         all_shared_lab_records = uploader.records_to_sync(apply_filters=False, full_sync=True).count()
    #         all_filtered_shared_lab_records = uploader.records_to_sync(apply_filters=True, full_sync=True).count()
    #         unsynced_filtered_shared_lab_records = uploader.records_to_sync(apply_filters=True, full_sync=False).count()
    #
    #         self.message_user(request, message=f"{destination} shared classifications for labs / passes filter / not yet synced : {all_shared_lab_records} / {all_filtered_shared_lab_records} / {unsynced_filtered_shared_lab_records}")

    @admin_action("Run now (delta sync)")
    def run_sync(self, request, queryset):
        self._run_sync(request, queryset)

    @admin_action("Run now (delta sync) single row")
    def run_sync_single(self, request, queryset):
        self._run_sync(request, queryset, max_rows=1)

    @admin_action("Run now (full sync)")
    def run_sync_full(self, request, queryset):
        sync_destination: SyncDestination
        for sync_destination in queryset:
            sync_destination.run(full_sync=True)
            self.message_user(request, message=f"Completed {str(sync_destination)}")


@admin.register(SyncRun)
class SyncRunAdmin(ModelAdminBasics):
    list_display = ('id', 'destination', 'created', 'status', 'full_sync', 'meta_short')
    list_filter = (('destination', RelatedFieldListFilter), 'status', 'full_sync')

    def has_add_permission(self, request):
        return False

    @admin_list_column(short_description="Meta", limit=100)
    def meta_short(self, obj: SyncRun):
        if meta := obj.meta:
            return json.dumps(meta)
        return ""

    @admin_action(short_description="Download report")
    def download_report(self, request, queryset: QuerySet[SyncRun]):
        if queryset.count() != 1:
            self.message_user(request, message="Can only run report on a single SyncRun at a time")
        else:
            sync_run: SyncRun = queryset.first()
            sync_runner = sync_runner_for_destination(sync_run.destination)
            return sync_runner.report_on(sync_run)


@admin.register(ClassificationModificationSyncRecord)
class ClassificationModificationSyncRecordAdmin(ModelAdminBasics):
    list_display = ('id', 'run', 'created', 'classification_modification', 'success', 'meta')
    search_fields = ('classification_modification__classification__lab_record_id', 'classification_modification__classification__id')

    def has_add_permission(self, request):
        return False
