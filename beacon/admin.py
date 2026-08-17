from django.contrib import admin

from beacon.models import BeaconInboundQuery, BeaconQueryCache


@admin.register(BeaconInboundQuery)
class BeaconInboundQueryAdmin(admin.ModelAdmin):
    list_display = ("created", "granularity", "authenticated", "observations_exists",
                    "observations_count", "classifications_exists", "classifications_count")
    list_filter = ("granularity", "authenticated", "observations_exists", "classifications_exists")
    date_hierarchy = "created"
    readonly_fields = ("created", "request_json")


@admin.register(BeaconQueryCache)
class BeaconQueryCacheAdmin(admin.ModelAdmin):
    list_display = ("variant", "node_id", "exists", "count", "error", "created")
    list_filter = ("node_id", "exists")
    search_fields = ("variant__id",)
    readonly_fields = ("created",)
