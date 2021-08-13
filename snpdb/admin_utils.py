import csv
import inspect
from typing import Optional

from django.contrib import admin, messages
from django.db import models
from django.db.models import AutoField, ForeignKey, DateTimeField
from django.http import HttpResponse
from django.utils.encoding import smart_text
# https://stackoverflow.com/questions/41228687/how-to-decorate-admin-actions-in-django-action-name-used-as-dict-key
from django_json_widget.widgets import JSONEditorWidget
from guardian.admin import GuardedModelAdmin


def admin_action(short_description: str):
    """
    Decorator, if used in ModelAdminBasics
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            # Just so the user knows that something happened no matter what
            request = args[1]
            queryset = args[2]
            messages.info(request, message=f"Called \"{short_description}\" on {queryset.count()} records")

            return method(*args, **kwargs)
        wrapper.short_description = short_description
        wrapper.__name__ = method.__name__
        wrapper.is_action = True
        return wrapper
    return decorator


def admin_list_column(short_description: Optional[str] = None, order_field: Optional[str] = None):

    def decorator(method):
        def wrapper(*args, **kwargs):
            # empty wrapper, we just want to modify short_description and mark as is_action
            return method(*args, **kwargs)
        if short_description:
            wrapper.short_description = short_description
        if order_field:
            wrapper.admin_order_field = order_field
        wrapper.__name__ = method.__name__
        wrapper.is_action = True
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

    def __new__(cls,  *args, **kwargs):
        actions = {}

        actions = dict((name, func) for name, func in inspect.getmembers(cls, lambda x: getattr(x, 'is_action', False)))
        if actions:
            cls.actions = ['export_as_csv'] + list(actions.keys())
        return super().__new__(cls)

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