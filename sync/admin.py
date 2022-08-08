from typing import Optional

from django.contrib import admin
from sync import models
from sync.models.models import SyncDestination
from sync.shariant.shariant_upload import ClassificationUploader


class ByDestinationFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'Destination Filter'
    parameter_name = 'destination'
    default_value = None

    def lookups(self, request, model_admin):
        return [(destination.pk, destination.name) for destination in SyncDestination.objects.all()]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(destination=self.value())
        return queryset


class SyncDestinationAdmin(admin.ModelAdmin):
    list_display = ('name', 'config', 'enabled')

    def _run_sync(self, request, queryset, max_rows: Optional[int] = None):
        sync_destination: SyncDestination
        for sync_destination in queryset:
            sync_destination.run(full_sync=False, max_rows=max_rows)
            self.message_user(request, message=f"Completed {str(sync_destination)} row limit = {max_rows}")

    def eligible_row_count(self, request, queryset):
        for destination in queryset:
            uploader = ClassificationUploader(destination)
            all_shared_lab_records = uploader.records_to_sync(apply_filters=False, full_sync=True).count()
            all_filtered_shared_lab_records = uploader.records_to_sync(apply_filters=True, full_sync=True).count()
            unsynced_filtered_shared_lab_records = uploader.records_to_sync(apply_filters=True, full_sync=False).count()

            self.message_user(request, message=f"{destination} shared classifications for labs / passes filter / not yet synced : {all_shared_lab_records} / {all_filtered_shared_lab_records} / {unsynced_filtered_shared_lab_records}")
    eligible_row_count.short_description = "Row Counts"

    def run_sync(self, request, queryset):
        self._run_sync(request, queryset)
    run_sync.short_description = "Run now (delta sync)"

    def run_sync_single(self, request, queryset):
        self._run_sync(request, queryset, max_rows=1)
    run_sync_single.short_description = "Run now (delta sync) single row"

    def run_sync_ten(self, request, queryset):
        self._run_sync(request, queryset, max_rows=10)
    run_sync_ten.short_description = "Run now (delta sync) 10 rows"

    def run_sync_full(self, request, queryset):
        sync_destination: SyncDestination
        for sync_destination in queryset:
            sync_destination.run(full_sync=True)
            self.message_user(request, message=f"Completed {str(sync_destination)}")
    run_sync_full.short_description = "Run now (full sync)"

    actions = [run_sync, run_sync_single, run_sync_ten, run_sync_full]


class SyncRunAdmin(admin.ModelAdmin):
    list_display = ('id', 'destination', 'created', 'status', 'meta')
    list_filter = (ByDestinationFilter,)


class ClassificationModificationSyncRecordAdmin(admin.ModelAdmin):
    list_display = ('id', 'run', 'created', 'classification_modification', 'success', 'meta')


admin.site.register(models.SyncDestination, SyncDestinationAdmin)
admin.site.register(models.SyncRun, SyncRunAdmin)
admin.site.register(models.ClassificationModificationSyncRecord, ClassificationModificationSyncRecordAdmin)
