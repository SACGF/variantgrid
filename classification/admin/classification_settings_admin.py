"""
Edits LabClassificationSettings / OrganizationClassificationSettings as inlines on snpdb's Lab and Organization admin:
the snpdb ModelAdmins are re-registered here with the inline added.
"""
from django.contrib import admin
from django.contrib.admin.widgets import AdminTextInputWidget
from django.db.models import TextField

from classification.models.classification_settings_models import (
    LabClassificationSettings,
    OrganizationClassificationSettings,
)
from snpdb.admin import LabAdmin, OrganizationAdmin
from snpdb.models.models import Lab, Organization


class LabClassificationSettingsInline(admin.StackedInline):
    model = LabClassificationSettings
    can_delete = False
    formfield_overrides = {
        TextField: {'widget': AdminTextInputWidget}
    }


class OrganizationClassificationSettingsInline(admin.StackedInline):
    model = OrganizationClassificationSettings
    can_delete = False


class LabClassificationAdmin(LabAdmin):
    list_display = (*LabAdmin.list_display, 'research', 'clinvar_key')
    inlines = (LabClassificationSettingsInline, )

    def get_queryset(self, request):
        return super().get_queryset(request).select_related('classification_settings__clinvar_key')

    @admin.display(boolean=True)
    def research(self, obj: Lab) -> bool:
        return obj.classification_settings.research

    def clinvar_key(self, obj: Lab):
        return obj.classification_settings.clinvar_key


class OrganizationClassificationAdmin(OrganizationAdmin):
    inlines = (OrganizationClassificationSettingsInline, )


admin.site.unregister(Lab)
admin.site.register(Lab, LabClassificationAdmin)
admin.site.unregister(Organization)
admin.site.register(Organization, OrganizationClassificationAdmin)
