from django.contrib import admin
from django.contrib.admin import RelatedFieldListFilter, TabularInline

from flags.models import Flag, FlagComment, FlagCollection
from flags.models.models import FlagType
from snpdb.admin_utils import ModelAdminBasics, AllValuesChoicesFieldListFilter


class FlagCommentAdminTabular(TabularInline):
    model = FlagComment

    fields = ("id", "created", "user", "text", "resolution")
    readonly_fields = ("created", "resolution",)

    def has_add_permission(self, request, obj):
        return False


@admin.register(Flag)
class FlagAdmin(ModelAdminBasics):
    search_fields = ('pk',)
    list_display = ('id', 'collection', 'flag_type', 'resolution', 'user', 'data', 'created', 'modified')
    list_filter = (('flag_type', RelatedFieldListFilter), ('resolution__status', AllValuesChoicesFieldListFilter), ('user', RelatedFieldListFilter))
    inlines = (FlagCommentAdminTabular,)

    def is_readonly_field(self, f) -> bool:
        if f.name == 'flag_type':
            return False
        return super().is_readonly_field(f)

    def has_add_permission(self, request):
        return False


@admin.register(FlagComment)
class FlagCommentAdmin(ModelAdminBasics):
    list_display = ('id', 'flag', 'user', 'text', 'resolution', 'created', 'modified')


@admin.register(FlagType)
class FlagTypeAdmin(ModelAdminBasics):
    list_display = ('id', 'context', 'label', 'description', 'help_text', 'raise_permission')
    list_filter = (('context', RelatedFieldListFilter), )


class FlagInline(TabularInline):
    model = Flag
    fields = ("id", "flag_type", "user", "resolution", "data")
    show_change_link = True

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
