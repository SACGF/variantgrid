from collections import defaultdict
import datetime
from typing import List, Tuple, Dict

from dateutil import parser
from django.conf import settings
from django.contrib import messages
from django.contrib.admin.views.decorators import staff_member_required
from django.core.exceptions import PermissionDenied, ValidationError, ObjectDoesNotExist
from django.db.models.aggregates import Count, Max
from django.db.models.base import ModelBase
from django.db.models.fields.reverse_related import OneToOneRel
from django.db.models.query_utils import Q
from django.urls.base import reverse_lazy
from django.utils import timezone
from django.utils.timezone import localtime
from functools import reduce
from functools import wraps, partial
import nameparser
import operator
from redis import Redis

from library.utils import invert_dict


def get_url_from_view_path(view_path):
    from django.contrib.sites.models import Site
    current_site = Site.objects.get_current()
    protocol = 'http'
    if 'shariant.org.au' in current_site.domain:
        protocol = 'https'
    return f'{protocol}://{current_site.domain}{view_path}'


def add_save_message(request, valid, name, created=False):
    action = "created" if created else "saved"

    if valid:
        msg = f"{name} {action} successfully."
        status = messages.INFO
    else:
        msg = f"Error. {name} not {action}."
        status = messages.ERROR

    messages.add_message(request, status, msg, extra_tags='save-message')


def require_superuser(f):

    @wraps(f)
    def wrapper(request, *args, **kwds):
        if request.user.is_superuser:
            return f(request, *args, **kwds)
        raise PermissionDenied("You must be a super user to view this page")

    return wrapper


def get_model_fields(model, ignore_fields=None):
    ignore_fields = set(ignore_fields or [])
    return [f.name for f in model._meta.fields if f.name not in ignore_fields]


def get_expanded_field(obj, field):
    """ Uses get_field_display if available """
    display_method = f"get_{field}_display"
    display_func = getattr(obj, display_method, None)
    if display_func:
        value = display_func()
    else:
        value = getattr(obj, field)
    return value


def get_model_fields_and_formatted_values_tuples_list(model):
    rows = []
    for name in get_model_fields(model):
        try:
            get_display_func = getattr(model, f"get_{name}_display")
            field = get_display_func()
        except:
            raw_field = getattr(model, name)
            if raw_field is not None:
#                logging.info("%s: %s", name, type(raw_field))
                if isinstance(raw_field, datetime.datetime):
                    raw_field = localtime(raw_field).strftime("%Y-%m-%d %H:%M:%S")

                field = str(raw_field)
            else:
                field = ''

        if field:
            rows.append((name, field))
    return rows


def column_arrays_from_values_queryset(qs, *fields, **formatters):
    column_arrays = defaultdict(list)
    keys = [f.split('__')[-1] for f in fields]
    annotation_run_values = qs.values(*fields)
    for value in annotation_run_values:
        for (k, field) in zip(keys, fields):
            v = value[field]
            f = formatters.get(field)
            if f:
                v = f(v)
            column_arrays[k].append(v)

    return column_arrays


def get_choices_formatter(choices, default=None):
    lookup = {}
    for choice, label in choices:
        lookup[choice] = label
        lookup[label] = label

    def get_choice(c):
        return lookup.get(c, default)

    return get_choice


def single_string_to_first_last_name_q(single_string):
    """ Returns a DJango Q object """
    name = nameparser.HumanName(single_string)

    queries = []
    first_name_queries = []
    first_name = name.first
    if first_name:
        first_name_queries.append(Q(first_name__iexact=first_name))

    last_name = name.last
    if last_name:
        queries.append(Q(last_name__iexact=last_name))
        # If we have the last name, also accept an initial
        first_name_queries.append(Q(first_name__istartswith=first_name[0]))

    queries.append(reduce(operator.or_, first_name_queries))
    return reduce(operator.and_, queries)


def count_values_qs(qs, field, count_field="count"):
    return qs.values(field).order_by(field).annotate(**{count_field: Count(field)})


def get_field_counts(qs, field):
    qs = count_values_qs(qs, field)
    field_counts = dict(qs.values_list(field, "count"))
    if None in field_counts:
        field_counts[None] = qs.filter(**{field + "__isnull": True}).count()
    return field_counts


staff_only = partial(staff_member_required, login_url=reverse_lazy('staff_only'))


def get_redis(**kwargs):
    port = kwargs.pop("port", settings.REDIS_PORT)
    return Redis(port=port, decode_responses=True, **kwargs)


def get_lower_choice(choices, value):
    d = invert_dict({name.lower(): v for (v, name) in choices})
    return d.get(value.lower())


def ensure_timezone_aware(datetime_date):
    aware_date = None
    if datetime_date:
        if isinstance(datetime_date, str):
            datetime_date = parser.parse(datetime_date)

        if timezone.is_aware(datetime_date):
            aware_date = datetime_date
        else:
            aware_date = timezone.make_aware(datetime_date, timezone.get_default_timezone())
    return aware_date


def ensure_mutally_exclusive_fields_not_set(obj, field_a, field_b):
    a_val = getattr(obj, field_a)
    b_val = getattr(obj, field_b)
    if a_val and b_val:
        msg = f"{obj} ({obj.pk}): You cannot set both fields '{field_a}' ('{a_val}') and '{field_b}' ('{b_val}') at the same time."
        raise ValueError(msg)


def set_form_read_only(form):
    for field in form.fields.values():
        field.disabled = True


def thread_safe_unique_together_get_or_create(klass, **kwargs):
    """ get_or_create has race conditions across threads.
        So, use unique_together to prevent a duplicate object being created, then
        fall back on normal get.

        If you get MultipleObjectsReturned - You need to add the unique_together
        constraint on your model to prevent creating duplicate objects """
    try:
        return klass.objects.get_or_create(**kwargs)
    except ValidationError:
        return klass.objects.get(**kwargs), False


def get_models_dict_by_column(klass, column='pk'):
    models_dict = {}
    for record in klass.objects.all():
        key = getattr(record, column, None)
        if key is not None:
            models_dict[key] = record
    return models_dict


def object_is_referenced(obj):
    return any(related_objects(obj))


def related_objects(obj):
    related_list = []

    one_to_one, many_to_many = object_relations(obj)
    for rel in one_to_one:
        try:
            related_object = getattr(obj, rel.get_accessor_name())
            related_list.append(related_object)
        except ObjectDoesNotExist:
            pass

    for rel in many_to_many:
        related_manager = getattr(obj, rel.get_accessor_name())
        if related_manager.exists():
            related_list.append(related_manager.all())

    return related_list


def object_relations(obj):
    """ returns (one_to_one, many_to_many) """

    is_1_to_1 = lambda field: isinstance(field, OneToOneRel)
    return discrimine(is_1_to_1, obj._meta._get_fields(forward=False))


def discrimine(pred, sequence):
    """Split a collection in two collections using a predicate.
    >>> discrimine(lambda x: x < 5, [3, 4, 5, 6, 7, 8])
    ... ([3, 4], [5, 6, 7, 8])
    """
    positive, negative = [], []
    for item in sequence:
        if pred(item):
            positive.append(item)
        else:
            negative.append(item)
    return positive, negative


class SortMetaOrderingMixin:
    """ Declares a '<' operator on Model - so you can sort lists same as querysets (driven by Meta.ordering) """
    def __lt__(self, other):
        for f in self._meta.ordering:
            v = getattr(self, f)
            ov = getattr(other, f)
            if v != ov:
                return v < ov
        return False


class SortByPKMixin:
    def __lt__(self, other):
        return self.pk < other.pk


def qs_highest_pk(qs):
    data = qs.aggregate(max_pk=Max("pk"))
    return data["max_pk"]


def highest_pk(model: ModelBase):
    qs = model.objects.all()
    return qs_highest_pk(qs)


def bulk_insert_class_data(apps, app_name: str, klass_name_and_data_list: List[Tuple[str, List[Dict]]]):
    """ For Django migrations """
    for klass_name, data in klass_name_and_data_list:
        klass = apps.get_model(app_name, klass_name)
        records = []
        for kwargs in data:
            records.append(klass(**kwargs))
        klass.objects.bulk_create(records)
