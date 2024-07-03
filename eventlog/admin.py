from django.contrib import admin
from django.contrib.auth.models import Group

from eventlog.models import ViewEvent
from snpdb.admin_utils import ModelAdminBasics
from snpdb.models import Lab, Organization


class AdminUsers(admin.SimpleListFilter):
    title = 'User Type'
    parameter_name = 'user_type'

    def lookups(self, request, model_admin):
        return (
            ('exclude', 'Exclude Admin/Test Users'),
        )

    def queryset(self, request, queryset):
        if self.value() == 'exclude':
            groups_to_exclude = ['variantgrid/tester', 'variantgrid/bot']
            groups = Group.objects.filter(name__in=groups_to_exclude)
            return queryset.exclude(user__groups__in=groups).exclude(user__is_superuser=True)
        return queryset


class LabFilter(admin.SimpleListFilter):
    title = 'Lab'
    parameter_name = 'lab'

    def lookups(self, request, model_admin):
        labs = Lab.objects.all().values_list('group_name', 'name')
        return [(group_name, name) for group_name, name in labs]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(user__groups__name=self.value())
        return queryset


class OrganizationFilter(admin.SimpleListFilter):
    title = 'Organization'
    parameter_name = 'organization'

    def lookups(self, request, model_admin):
        organizations = Organization.objects.all().values_list('group_name', 'name')
        return [(group_name, name) for group_name, name in organizations]

    def queryset(self, request, queryset):
        if self.value():
            return queryset.filter(user__groups__name=self.value())
        return queryset


@admin.register(ViewEvent)
class ViewEventAdmin(ModelAdminBasics):

    """
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    view_name = models.TextField()
    args = models.JSONField(null=False, blank=True, default=dict)
    path = models.TextField()
    method = models.TextField()

    """

    list_display = ('pk', 'view_name', 'user', 'args', 'path', 'created')
    search_fields = ('view_name', 'path')
    list_filter = (
        OrganizationFilter,
        LabFilter,
        AdminUsers
    )
