from typing import Optional

from django.contrib import admin
from django.contrib.admin import RelatedFieldListFilter

from snpdb.admin_utils import ModelAdminBasics, admin_action
from sync.models import SyncRun, ClassificationModificationSyncRecord
from sync.models.models import SyncDestination


@admin.register(SyncDestination)
class SyncDestinationAdmin(ModelAdminBasics):
    list_display = ('name', 'config', 'enabled')

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget()
        }, **kwargs)

    def _run_sync(self, request, queryset, max_rows: Optional[int] = None):
        sync_destination: SyncDestination
        for sync_destination in queryset:
            sync_destination.run(full_sync=False, max_rows=max_rows)
            self.message_user(request, message=f"Completed {str(sync_destination)} row limit = {max_rows}")

    @admin_action("Run now (delta sync)")
    def run_sync(self, request, queryset):
        self._run_sync(request, queryset)

    @admin_action("Run now (delta sync) single row")
    def run_sync_single(self, request, queryset):
        self._run_sync(request, queryset, max_rows=1)

    @admin_action("Run now (delta sync) 10 rows")
    def run_sync_ten(self, request, queryset):
        self._run_sync(request, queryset, max_rows=10)

    @admin_action("Run now (full sync)")
    def run_sync_full(self, request, queryset):
        sync_destination: SyncDestination
        for sync_destination in queryset:
            sync_destination.run(full_sync=True)
            self.message_user(request, message=f"Completed {str(sync_destination)}")


@admin.register(SyncRun)
class SyncRunAdmin(ModelAdminBasics):
    list_display = ('id', 'destination', 'created', 'status', 'meta')
    list_filter = (('destination', RelatedFieldListFilter),)

    def has_add_permission(self, request):
        return False


@admin.register(ClassificationModificationSyncRecord)
class ClassificationModificationSyncRecordAdmin(ModelAdminBasics):
    list_display = ('id', 'run', 'created', 'classification_modification', 'success', 'meta')
    search_fields = ('classification_modification__classification__lab_record_id', 'classification_modification__classification__id')

    def has_add_permission(self, request):
        return False
