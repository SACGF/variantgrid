from django.contrib import admin
from django.contrib.admin import RelatedFieldListFilter, TabularInline

from flags.models import Flag, FlagComment, FlagTypeResolution, FlagCollection
from flags.models.models import FlagType
from snpdb.admin_utils import ModelAdminBasics, AllValuesChoicesFieldListFilter


class FlagCommentAdmin(TabularInline):
    model = FlagComment

    def has_add_permission(self, request, obj):
        return False


@admin.register(Flag)
class FlagAdmin(ModelAdminBasics):
    list_display = ('id', 'collection', 'flag_type', 'resolution', 'user', 'data', 'created', 'modified')
    list_filter = (('flag_type', RelatedFieldListFilter), ('resolution__status', AllValuesChoicesFieldListFilter), ('user', RelatedFieldListFilter))
    inlines = (FlagCommentAdmin,)

    def has_add_permission(self, request):
        return False


class FlagTypeResolution(TabularInline):
    model = FlagTypeResolution

    def has_add_permission(self, request, obj):
        return False

    def has_delete_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False

@admin.register(FlagType)
class FlagTypeAdmin(ModelAdminBasics):
    list_display = ('id', 'context', 'label', 'description', 'help_text', 'raise_permission')
    list_filter = (('context', RelatedFieldListFilter), )
    inlines = (FlagTypeResolution,)


class FlagInline(TabularInline):
    model = Flag

    def has_add_permission(self, request, obj):
        return False

    def has_delete_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False


@admin.register(FlagCollection)
class FlagCollectionAdmin(ModelAdminBasics):
    list_display = ('pk', 'text')
    list_filter = (('context', RelatedFieldListFilter), )
    search_fields = ('pk',)
    inlines = (FlagInline,)

    def text(self, obj: FlagCollection):
        return str(obj)

    def has_add_permission(self, request):
        return False
