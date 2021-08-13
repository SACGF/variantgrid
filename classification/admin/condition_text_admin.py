from django.contrib import admin
from django.contrib.admin import RelatedFieldListFilter
from django.db.models import QuerySet
from classification.models import ConditionText, ConditionTextMatch
from library.guardian_utils import admin_bot
from snpdb.admin_utils import ModelAdminBasics, admin_action


class ConditionTextStatusFilter(admin.SimpleListFilter):
    title = 'Status Filter'
    parameter_name = 'status_filter'
    default_value = None

    def lookups(self, request, model_admin):
        return [("outstanding", "Only with outstanding classifications")]

    def queryset(self, request, queryset):
        if value := self.value():
            if value == "outstanding":
                queryset = queryset.exclude(classifications_count_outstanding=0)
        return queryset


class ConditionTextMatchAdmin(admin.TabularInline):
    model = ConditionTextMatch

    def has_add_permission(self, request, obj=None):
        return False

    def has_change_permission(self, request, obj=None):
        return False


@admin.register(ConditionText)
class ConditionTextAdmin(ModelAdminBasics):
    inlines = (ConditionTextMatchAdmin,)
    search_fields = ('id', 'normalized_text')
    list_display = ["pk", "lab", "normalized_text", "classifications_count", "classifications_count_outstanding"]
    list_filter = (ConditionTextStatusFilter, ('lab', RelatedFieldListFilter))

    def has_add_permission(self, request):
        return False

    def get_form(self, request, obj=None, **kwargs):
        return super().get_form(request, obj, widgets={
            'normalized_text': admin.widgets.AdminTextInputWidget()
        }, **kwargs)

    @admin_action("Automatch (leaves existing data)")
    def auto_match(self, request, queryset: QuerySet[ConditionText]):
        for condition_text in queryset:
            ConditionTextMatch.attempt_automatch(condition_text=condition_text)
            condition_text.save()

    @admin_action("Clear any matches")
    def clear(self, request, queryset):
        condition_text: ConditionText
        for condition_text in queryset:
            condition_text.clear()


class ConditionTextMatchUserFilter(admin.SimpleListFilter):
    list_per_page = 200
    title = 'User Filter'
    parameter_name = 'user'
    default_value = None

    def lookups(self, request, model_admin):
        return [
            ("any", "Any User"),
            ("non_admin", "Non Admin"),
            ("bot", "Auto-Assigned")
        ]

    def queryset(self, request, queryset):
        if user := self.value():
            if user == "any":
                queryset = queryset.filter(last_edited_by__isnull=False).exclude(last_edited_by=admin_bot())
            elif user == "non_admin":
                queryset = queryset.filter(last_edited_by__is_superuser=False).exclude(last_edited_by=admin_bot())
            elif user == "bot":
                queryset = queryset.filter(last_edited_by=admin_bot())
        return queryset
