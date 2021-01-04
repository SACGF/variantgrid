from django.contrib import admin
from unidecode import unidecode
from snpdb.models import Organization, Lab
import re


def make_code_friendly(text: str) -> str:
    """
    convert accented characters to non-accented counterparts
    lower case, replace - and spaces with underscores
    remove anything that's then not a-z or underscore
    """
    text = unidecode(text) \
        .lower() \
        .replace('-', '_').replace(' ', '_')
    return re.sub(r'[^a-z0-9_]', '', text)


class LabAdmin(admin.ModelAdmin):
    list_per_page = 200
    list_display = ('name', 'group_name', 'organization', 'external', 'classification_config')

    fieldsets = (
        ('Basic', {'fields': ('name', 'group_name', 'organization')}),
        ('Position', {'fields': ('city', 'state', 'country', 'lat', 'long')}),
        ('Style', {'fields': ('url', 'css_class')}),
        ('Submissions', {'fields': ('classification_config', 'upload_location', 'external')})
    )

    def get_form(self, request, obj=None, **kwargs):

        return super(LabAdmin, self).get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'institution': admin.widgets.AdminTextInputWidget(),
            'group_name': admin.widgets.AdminTextInputWidget(),
            'city': admin.widgets.AdminTextInputWidget(),
            'state': admin.widgets.AdminTextInputWidget(),
            'country': admin.widgets.AdminTextInputWidget(),
            'lat': admin.widgets.AdminTextInputWidget(),
            'long': admin.widgets.AdminTextInputWidget(),
            'url': admin.widgets.AdminURLFieldWidget(),
            'css_class': admin.widgets.AdminTextInputWidget(),
            'upload_location': admin.widgets.AdminTextInputWidget()

        }, **kwargs)

    def fix_group_name(self, request, queryset):
        safety_reg = re.compile(r'^[a-z0-9_]*$')
        fixed = 0
        already_good = 0
        lab: Lab
        for lab in queryset:
            org_group_name = lab.organization.group_name
            if not lab.group_name or not safety_reg.match(lab.group_name) or not lab.group_name.startswith(org_group_name):
                lab.group_name = org_group_name + '/' + make_code_friendly(lab.name)
                lab.save()
                fixed = fixed + 1
            else:
                already_good = already_good + 1
        self.message_user(request, f"{fixed} updated, {already_good} no change required")

    fix_group_name.short_description = 'Fix group name'

    actions = [fix_group_name]


class OrganizationAdmin(admin.ModelAdmin):

    list_display = ('name', 'group_name', 'classification_config')

    fieldsets = (
        ('Basic', {'fields': ('name', 'short_name', 'group_name', 'active')}),
        ('Submissions', {'fields': ('classification_config', 'classification_report_template')})
    )

    def fix_group_name(self, request, queryset):
        org: Organization
        safety_reg = re.compile(r'^[a-z0-9_]*$')
        fixed = 0
        already_good = 0
        for org in queryset:
            if not org.group_name or not safety_reg.match(org.group_name):
                org.group_name = make_code_friendly(org.name)
                org.save()
                fixed = fixed + 1
            else:
                already_good = already_good + 1
        self.message_user(request, f"{fixed} updated, {already_good} no change required")

    fix_group_name.short_description = 'Fix group name'

    actions = [fix_group_name]

    def get_form(self, request, obj=None, **kwargs):
        return super(OrganizationAdmin, self).get_form(request, obj, widgets={
            'name': admin.widgets.AdminTextInputWidget(),
            'short_name': admin.widgets.AdminTextInputWidget(),
            'group_name': admin.widgets.AdminTextInputWidget()

        }, **kwargs)

# The JSONEditor has a bug in it that stops patternProperties
# from being useful. Specifically if you try to add arbitrary properties
# you get locked in the UI and can't continue to edit.

# classification_config_schema_preferred = {
#     'type': 'object',
#         'patternProperties': {
#         "[a-zA-Z0-9_]+": {
#             "oneOf": [
#                 { "type": "boolean" },
#                 {
#                     "type": 'object',
#                     'properties': {
#                         'hide': {'type': 'boolean' },
#                         'custom_options': {
#                             'type': 'array', 'items': {
#                                 'type': 'string'
#                             }
#                         }
#                     }
#                 }
#             ]
#         }
#     }
# }
