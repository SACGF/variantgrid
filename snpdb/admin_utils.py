import csv

from django.contrib import admin
from django.db import models
from django.db.models import AutoField, ForeignKey, DateTimeField
from django.http import HttpResponse
from django.utils.encoding import smart_text
# https://stackoverflow.com/questions/41228687/how-to-decorate-admin-actions-in-django-action-name-used-as-dict-key
from django_json_widget.widgets import JSONEditorWidget
from guardian.admin import GuardedModelAdmin


def short_description(short_description: str):
    def decorator(admin_action):
        def wrapper(*args, **kwargs):
            return admin_action(*args, **kwargs)
        wrapper.short_description = short_description
        wrapper.__name__ = short_description
        return wrapper
    return decorator


class AllValuesChoicesFieldListFilter(admin.AllValuesFieldListFilter):

    def choices(self, changelist):
        yield {
            'selected': self.lookup_val is None and self.lookup_val_isnull is None,
            'query_string': changelist.get_query_string({}, [self.lookup_kwarg, self.lookup_kwarg_isnull]),
            'display': 'All',
        }
        include_none = False

        # all choices for this field
        choices = dict(self.field.choices)

        for val in self.lookup_choices:
            if val is None:
                include_none = True
                continue
            val = smart_text(val)
            yield {
                'selected': self.lookup_val == val,
                'query_string': changelist.get_query_string({
                    self.lookup_kwarg: val,
                }, [self.lookup_kwarg_isnull]),

                # instead code, display title
                'display': choices[val],
            }
        if include_none:
            yield {
                'selected': bool(self.lookup_val_isnull),
                'query_string': changelist.get_query_string({
                    self.lookup_kwarg_isnull: 'True',
                }, [self.lookup_kwarg]),
                'display': self.empty_value_display,
            }


class AdminExportCsvMixin:
    def export_as_csv(self, request, queryset):

        meta = self.model._meta
        field_names = [field.name for field in meta.fields]

        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={}.csv'.format(meta)
        writer = csv.writer(response)

        writer.writerow(field_names)
        for obj in queryset:
            writer.writerow([getattr(obj, field) for field in field_names])
        return response

    export_as_csv.short_description = "Export Selected as CSV"

    def _is_readonly(self, f) -> bool:
        if not f.editable:
            return True  # does this make all the below redundant?
        if isinstance(f, (AutoField, ForeignKey)):
            return True
        if isinstance(f, DateTimeField):
            if f.auto_now or f.auto_now_add:
                return True
        return False

    def _get_readonly_fields(self, request, obj=None):
        return [f.name for f in self.model._meta.fields if self._is_readonly(f)]

    def _get_fields(self, request, obj=None, **kwargs):
        first = []
        second = []
        for f in self.model._meta.fields:
            if isinstance(f, (AutoField, ForeignKey)):
                # put ids and foreign keys first
                # first.append(f.name)
                first.append(f.name)
            else:
                second.append(f.name)

        return first + second


class ModelAdminBasics(admin.ModelAdmin, AdminExportCsvMixin):
    formfield_overrides = {
        models.JSONField: {'widget': JSONEditorWidget},
    }

    # wanted to call this BaseModelAdmin but that was already taken
    actions = ["export_as_csv"]

    def get_readonly_fields(self, request, obj=None):
        return self._get_readonly_fields(request=request, obj=obj)

    def get_fields(self, request, obj=None):
        return self._get_fields(request=request, obj=obj)


class GuardedModelAdminBasics(GuardedModelAdmin, AdminExportCsvMixin):
    actions = ["export_as_csv"]

    def get_readonly_fields(self, request, obj=None):
        return self._get_readonly_fields(request=request, obj=obj)

    def get_fields(self, request, obj=None):
        return self._get_fields(request=request, obj=obj)