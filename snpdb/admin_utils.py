import inspect
from functools import cached_property
from typing import Optional, Iterator, Type

from dateutil.tz import gettz
from django.conf import settings
from django.contrib import admin, messages
from django.db import models
from django.db.models import AutoField, ForeignKey, DateTimeField, Model
from django.forms import Widget
from django.http import StreamingHttpResponse, HttpResponseRedirect
from django.http.response import HttpResponseBase
from django.urls import path, NoReverseMatch, reverse
from django.utils.encoding import smart_str
from django.utils.safestring import SafeString
from django_json_widget.widgets import JSONEditorWidget
from guardian.admin import GuardedModelAdminMixin

from library.log_utils import log_admin_change
from library.utils import delimited_row, WrappablePartial, limit_str


class AllValuesChoicesFieldListFilter(admin.AllValuesFieldListFilter):
    """
    Used to provide filter for choice fields (will display the choice label and not just the code)
    Unlike the default filter for a regular field, it will provide all TextChoices and not just
    all distinct values
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


def admin_action(short_description: str):
    """
    Decorator, if used in ModelAdminBasics, marks function as something that can be done on selected rows.
    As a bonus will also send an info message confirming to the user how many rows were selected for the action
    From https://stackoverflow.com/questions/41228687/how-to-decorate-admin-actions-in-django-action-name-used-as-dict-key
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

            for record in queryset:
                log_admin_change(record, short_description)
            return method(*args, **kwargs)
        wrapper.short_description = short_description
        wrapper.line_number = inspect.getsourcelines(method)[1]
        wrapper.__name__ = method.__name__
        wrapper.is_action = True
        return wrapper
    return decorator


def admin_model_action(url_slug: str, short_description: Optional[str] = None, icon: Optional[str] = None):
    """
    Decorator, if used in ModelAdminBasics, marks function as top level across the model, rather than for any specific record
    """
    def decorator(method):
        def wrapper(*args, **kwargs):
            # empty wrapper, we just want to modify short_description and mark as is_action
            result = method(*args, **kwargs)
            if not result:
                result = HttpResponseRedirect("../")
            return result

        wrapper.method = method
        wrapper.url_slug = url_slug
        wrapper.short_description = short_description
        wrapper.line_number = inspect.getsourcelines(method)[1]
        wrapper.icon = icon
        wrapper.__name__ = method.__name__
        wrapper.is_model_action = True
        return wrapper

    return decorator


def admin_list_column(
        short_description: Optional[str] = None,
        order_field: Optional[str] = None,
        is_boolean: bool = False,
        limit: int = 100
        ):
    """
    Decorator, mark a function as acting like an admin list column if class extends ModelAdminBasics
    Note that the column still needs to be added to "list_display" because its' hard to order the data otherwise
    :short_description: Optional Overrides the default display name
    :order_field: Optional, If provided is used to order the table using this database column
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            # empty wrapper, we just want to modify short_description and mark as is_action
            result = method(*args, **kwargs)
            if result and limit and not isinstance(result, (bool, int, float)):
                if not isinstance(result, SafeString):
                    # don't truncate SafeStrings as would include HTML that could easily break
                    result = limit_str(str(result), limit)
            return result

        if short_description:
            wrapper.short_description = short_description
        if order_field:
            wrapper.admin_order_field = order_field
        if is_boolean:
            wrapper.boolean = True
        wrapper.__name__ = method.__name__
        wrapper.is_list_column = True
        return wrapper
    return decorator


def inject_self_to_decorated(wrapped, self_instance):
    """
    Reflectively called the wrapped method, injecting self, if this first argument of the method is self
    """
    takes_self = False
    if method_args := inspect.getfullargspec(wrapped.method).args:
        takes_self = method_args[0] == 'self'
    if takes_self:
        return WrappablePartial(wrapped, self_instance)
    else:
        return wrapped


def export_as_csv(modeladmin, request, queryset) -> HttpResponseBase:
    """
    Action to be provided against all models
    """
    meta = queryset.model._meta
    field_names = [field.name for field in meta.fields]
    if related_fields := [field.name for field in meta.fields if isinstance(field, ForeignKey)]:
        queryset = queryset.select_related(*related_fields)

    def get_formatted_attr(qs_obj2, field):
        display_name = f"get_{field}_display"
        if hasattr(qs_obj2, display_name):
            return getattr(qs_obj2, display_name)()
        else:
            return getattr(qs_obj2, field)

    def data_generator() -> Iterator[str]:
        nonlocal field_names
        nonlocal queryset

        yield delimited_row(field_names)

        for qs_obj in queryset:
            yield delimited_row([get_formatted_attr(qs_obj, field) for field in field_names])

    response = StreamingHttpResponse(data_generator(), content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename={meta}.csv'
    return response


export_as_csv.short_description = "Export selected as CSV"

# Sets this action to run on every single admin screen
admin.site.add_action(export_as_csv, 'export_as_csv')


class ModelAdminBasics(admin.ModelAdmin):
    """
    Make every Admin class extend this (or GuardedModelAdminBasics)
    *
    * Provides support for annotating methods @admin_method
    * Comes with export_as_csv for selected rows
    * Sets a nice default editor for JSONFields
    * Marks foreign keys as readonly - to prevent generating a million record select
        (override is_readonly_field and set value for autocomplete_fields or raw_id_fields
        if you need editable foreign fields)
    """

    formfield_overrides = {
        models.JSONField: {'widget': JSONEditorWidget},
        # sadly this doesn't work as general widget since "rel" and "admin_site" need values
        # models.ForeignKey: {'widget': ForeignKeyRawIdWidget}
    }

    def __new__(cls, *args, **kwargs):
        """
        Sets up actions (methods annotated with @admin_action)
        """
        actions = [func for _, func in inspect.getmembers(cls, lambda x: getattr(x, 'is_action', False))]
        if actions:
            actions.sort(key=lambda x: x.line_number)
        cls.actions = actions

        instance = super().__new__(cls)

        model_actions = [func for _, func in inspect.getmembers(cls, lambda x: getattr(x, 'is_model_action', False))]
        cls.model_urls = []
        cls.model_actions = model_actions
        if model_actions:
            model_actions.sort(key=lambda x: x.line_number)
            cls.model_urls = [path(ma.url_slug, inject_self_to_decorated(ma, instance)) for ma in model_actions]

        return instance

    def changelist_view(self, request, extra_context=None):
        # provides data to the template so that model wide methods can be added to the top
        extra_context = extra_context or {}
        extra_context["model_actions"] = [{
            "label": ma.short_description,
            "url": ma.url_slug,
            "icon": ma.icon or "fa-solid fa-play"
        } for ma in self.model_actions]
        return super().changelist_view(request, extra_context=extra_context)

    def changeform_view(self, request, object_id=None, form_url="", extra_context=None):
        # provides data to the template so that actions can be added to the top
        extra_context = extra_context or {}
        all_actions = self.get_actions(request)
        extra_context["actions"] = [{
            "id": a[1],
            "label": a[2]
        } for a in all_actions.values() if not a[1] == 'delete_selected']
        return super().changeform_view(request, object_id, form_url, extra_context)

    def single_page_action(self, request, object_id):
        # called by the change template when choosing a single action
        action_name = request.POST.get("action")
        if isinstance(object_id, str):
            # really tried to find the proper way to decode the object_id
            # but just left doing this sub, I'm sure this will fail for other IDs with harder to escape values
            object_id = object_id.replace("_5F", "_")
            object_id = object_id.replace("_3A", ":")

        action = self.get_action(action_name)
        response = action[0](self, request, self.model.objects.filter(pk=object_id))
        if not response:
            response = HttpResponseRedirect("../change/")
        return response

    def get_urls(self):
        # injects single_action, model wide URLs and thent he default URL generation
        resulting_urls = [path('<path:object_id>/single_action/', self.single_page_action)] + self.model_urls + super().get_urls()
        return resulting_urls

    @cached_property
    def tz(self):
        return gettz(settings.TIME_ZONE)

    def format_datetime(self, datetime) -> str:
        default_timezoned = datetime.astimezone(self.tz)
        return f"{default_timezoned.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3]}"

    def format_datetime_seconds(self, datetime) -> str:
        default_timezoned = datetime.astimezone(self.tz)
        return default_timezoned.strftime('%Y-%m-%d %H:%M:%S')

    def is_readonly_field(self, f) -> bool:
        # override to make fields editable
        if not f.editable:
            return True
        if isinstance(f, (AutoField, ForeignKey)):
            if f.name in self.autocomplete_fields:
                return False
            return True
        if isinstance(f, DateTimeField):
            if f.auto_now or f.auto_now_add:
                return True
        return False

    def _get_readonly_fields(self, request, obj=None) -> list[str]:
        return [f.name for f in self.model._meta.fields if self.is_readonly_field(f)]

    def _get_fields(self, request, obj=None, **kwargs) -> list[str]:
        first: list[str] = []
        second: list[str] = []
        exclude_us = self.exclude or set()
        for f in self.model._meta.fields:
            if f.name in exclude_us:
                continue
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

    def widget_overrides(self) -> dict[str, Widget]:
        return {}

    def get_form(self, request, obj=None, **kwargs):
        kwarg_widgets = kwargs.pop('widgets', {})
        overrides = self.widget_overrides()
        widgets = {}
        widgets.update(kwarg_widgets)
        widgets.update(overrides)

        return super().get_form(request, obj, widgets=widgets, **kwargs)


class GuardedModelAdminBasics(GuardedModelAdminMixin, ModelAdminBasics):
    pass


def get_admin_url(obj: Model):
    try:
        model = obj._meta.model
        meta = model._meta
        path = f"admin:{meta.app_label}_{meta.model_name}_change"
        return reverse(path, kwargs={"object_id": obj.pk})
    except NoReverseMatch:
        return None

def get_admin_model_url(model_type: Type[Model]):
    try:
        meta = model_type._meta
        path = f"admin:{meta.app_label}_{meta.model_name}_changelist"
        return reverse(path)
    except NoReverseMatch:
        return None
