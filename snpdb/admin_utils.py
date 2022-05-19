import inspect
from typing import Optional, List, Iterator

import pytz
from django.conf import settings
from django.contrib import admin, messages
from django.db import models
from django.db.models import AutoField, ForeignKey, DateTimeField
from django.http import HttpResponse, StreamingHttpResponse
from django.utils.encoding import smart_str
from django_json_widget.widgets import JSONEditorWidget
from guardian.admin import GuardedModelAdminMixin

# https://stackoverflow.com/questions/41228687/how-to-decorate-admin-actions-in-django-action-name-used-as-dict-key
from lazy import lazy

from library.utils import delimited_row


def admin_action(short_description: str):
    """
    Decorator, if used in ModelAdminBasics, marks function as something that can be done on selected rows.
    As a bonus will also send an info message confirming to the user how many rows were selected for the action
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            # Just so the user knows that something happened no matter what
            request = args[1]
            queryset = args[2]
            count = queryset.count()
            if count == 1:
                messages.info(request, message=f"Called \"{short_description}\" on \"{queryset.first()}\"")
            else:
                messages.info(request, message=f"Called \"{short_description}\" on {queryset.count()} records")

            return method(*args, **kwargs)
        wrapper.short_description = short_description
        wrapper.line_number = inspect.getsourcelines(method)[1]
        wrapper.__name__ = method.__name__
        wrapper.is_action = True
        return wrapper
    return decorator


def admin_list_column(short_description: Optional[str] = None, order_field: Optional[str] = None):
    """
    Mark a function as acting like an admin list column
    Note that the column still needs to be added to "list_display" as that is the easiest
    way to allow ordering of columns.
    :short_description: Optional Overrides the default display name
    :order_field: Optional, If provided is used to order the table using this database column
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            # empty wrapper, we just want to modify short_description and mark as is_action
            return method(*args, **kwargs)
        if short_description:
            wrapper.short_description = short_description
        if order_field:
            wrapper.admin_order_field = order_field
        wrapper.__name__ = method.__name__
        wrapper.is_list_column = True
        return wrapper
    return decorator


class AllValuesChoicesFieldListFilter(admin.AllValuesFieldListFilter):
    """
    Used to provide filter for choice fields (will display the choice label and not just the code)
    e.g. list_filter = ('status', AllValuesChoicesFieldListFilter)
    """

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
            val = smart_str(val)
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


class ModelAdminBasics(admin.ModelAdmin):
    """
    Make every Admin class extend this (or GuardedModelAdminBasics)
    *
    * Provides support for annotating methods @admin_method
    * Comes with export_as_csv for selected rows
    * Sets a nice default editor for JSONFields
    * Marks foreign keys as readonly (override is_readonly_field and set value for
    * autocomplete_fields or raw_id_fields if you need editable foreign fields)
    """

    formfield_overrides = {
        models.JSONField: {'widget': JSONEditorWidget},
        # sadly this doesn't work as general widget since "rel" and "admin_site" need values
        # models.ForeignKey: {'widget': ForeignKeyRawIdWidget}
    }

    def export_as_csv(self, request, queryset) -> HttpResponse:
        meta = self.model._meta
        field_names = [field.name for field in meta.fields]
        if related_fields := [field.name for field in meta.fields if isinstance(field, ForeignKey)]:
            queryset = queryset.select_related(*related_fields)

        def data_generator() -> Iterator[str]:
            nonlocal field_names
            nonlocal queryset

            yield delimited_row(field_names)
            for qs_obj in queryset:
                yield delimited_row([getattr(qs_obj, field) for field in field_names])

        response = StreamingHttpResponse(data_generator(), content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename={}.csv'.format(meta)
        return response

    export_as_csv.short_description = "Export selected as CSV"

    def __new__(cls, *args, **kwargs):
        """
        Sets up actions (methods annotated with @admin_action)
        """
        actions = [func for _, func in inspect.getmembers(cls, lambda x: getattr(x, 'is_action', False))]

        if actions:
            actions.sort(key=lambda x: x.line_number)

        cls.actions = ['export_as_csv'] + actions
        return super().__new__(cls)

    @lazy
    def tz(self):
        return pytz.timezone(settings.TIME_ZONE)

    def format_datetime(self, datetime) -> str:
        default_timezoned = datetime.astimezone(self.tz)
        return f"{default_timezoned.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]}"

    def is_readonly_field(self, f) -> bool:
        if not f.editable:
            return True  # does this make all the below redundant?
        if isinstance(f, (AutoField, ForeignKey)):
            return True
        if isinstance(f, DateTimeField):
            if f.auto_now or f.auto_now_add:
                return True
        return False

    def _get_readonly_fields(self, request, obj=None) -> List[str]:
        return [f.name for f in self.model._meta.fields if self.is_readonly_field(f)]

    def _get_fields(self, request, obj=None, **kwargs) -> List[str]:
        first: List[str] = list()
        second: List[str] = list()
        for f in self.model._meta.fields:
            if isinstance(f, (AutoField, ForeignKey)):
                # put ids and foreign keys first
                # first.append(f.name)
                first.append(f.name)
            else:
                second.append(f.name)

        return first + second

    def get_readonly_fields(self, request, obj=None):
        return self._get_readonly_fields(request=request, obj=obj)

    def get_fields(self, request, obj=None):
        return self._get_fields(request=request, obj=obj)


class GuardedModelAdminBasics(GuardedModelAdminMixin, ModelAdminBasics):
    pass
