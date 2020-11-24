from django.contrib import admin
from flags import models
from flags.models import FlagTypeResolution, FlagStatus
from flags.models.models import FlagType

class FlagTypeFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'Flag Type'
    parameter_name = 'flag_type'
    default_value = None

    def lookups(self, request, model_admin):
        return [(flag_type.id, flag_type.label) for flag_type in FlagType.objects.all()]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(flag_type=self.value())
        return queryset


class FlagStatusFilter(admin.SimpleListFilter):
    list_per_page = 50
    title = 'Resolution'
    parameter_name = 'resolution'
    default_value = None

    def lookups(self, request, model_admin):
        return FlagStatus.CHOICES

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(resolution__status=self.value())
        return queryset


class FlagAdmin(admin.ModelAdmin):
    list_display = ('id', 'collection', 'flag_type', 'resolution', 'user', 'data', 'created', 'modified')
    list_filter = (FlagTypeFilter, FlagStatusFilter,)

admin.site.register(models.FlagType)
admin.site.register(models.FlagTypeContext)
admin.site.register(models.FlagCollection)
admin.site.register(models.Flag, FlagAdmin)
admin.site.register(models.FlagResolution)
admin.site.register(models.FlagTypeResolution)
