from django.contrib import admin
from sync import models
from sync.models.models import SyncDestination


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
    list_display = ('name', 'config')

    def run_sync(self, request, queryset):
        for sync_destination in queryset:
            sync_destination.run(full_sync=False)
            self.message_user(request, message=f"Completed {str(sync_destination)}")
    run_sync.short_description = "Run now (delta sync)"

    def run_sync_full(self, request, queryset):
        for sync_destination in queryset:
            sync_destination.run(full_sync=True)
            self.message_user(request, message=f"Completed {str(sync_destination)}")
    run_sync_full.short_description = "Run now (full sync)"

    actions = [run_sync, run_sync_full]


class SyncRunAdmin(admin.ModelAdmin):
    list_display = ('id', 'destination', 'created')
    list_filter = (ByDestinationFilter,)


class ClassificationModificationSyncRecordAdmin(admin.ModelAdmin):
    list_display = ('id', 'run', 'created', 'classification_modification')


admin.site.register(models.SyncDestination, SyncDestinationAdmin)
admin.site.register(models.SyncRun, SyncRunAdmin)
admin.site.register(models.ClassificationModificationSyncRecord, ClassificationModificationSyncRecordAdmin)
