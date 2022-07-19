import csv
import difflib
import hashlib
import importlib
import inspect
import io
import json
import logging
import math
import operator
import re
import string
import subprocess
import time
import uuid
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta
from decimal import Decimal
from enum import Enum
from functools import reduce
from itertools import islice
from json.encoder import JSONEncoder
from operator import attrgetter
from typing import TypeVar, Optional, Iterator, Tuple, Any, List, Iterable, Set, Dict, Union, Callable, Type, Generic
from urllib.parse import urlparse
from django.utils.timezone import localtime
from bs4 import BeautifulSoup
from dateutil import parser
from django.conf import settings
from django.core.serializers import serialize
from django.db import models
from django.db.models.query import QuerySet
from django.http import StreamingHttpResponse, HttpRequest
from django.utils import html, timezone
from django.utils.functional import SimpleLazyObject
from django.utils.safestring import SafeString, mark_safe

from uicore.json.json_types import JsonObjType

FLOAT_REGEX = r'([-+]?[0-9]*\.?[0-9]+.|Infinity)'


# isn't restrictive enough, but helps documention
JSON = Union[Dict[str, Any], List[Any], int, float, str, bool]


class DjangoJSONEncoder(JSONEncoder):

    def default(self, obj):
        if isinstance(obj, QuerySet):
            # `default` must return a python serializable
            # structure, the easiest way is to load the JSON
            # string produced by `serialize` and return it
            return json.loads(serialize('json', obj))
        return JSONEncoder.default(self, obj)


def string_deterministic_hash(s: str) -> int:
    """
    When you want the same string to always hash to the same value (where security
    isnt a concern). Hashlib seemed overkill for this purpose.
    """
    val = 0
    for c in s:
        val = val * 13 + ord(c)
    return val


def md5sum_str(s):
    s_bytes = s.encode()
    return hashlib.md5(s_bytes).hexdigest()


def sha1_str(s):
    s_bytes = s.encode()
    return hashlib.sha1(s_bytes).hexdigest()


def invert_dict(source_dict):
    return {v: k for k, v in source_dict.items()}


def invert_dict_of_lists(source_dict):
    d = {}
    for k, v_list in source_dict.items():
        for v in v_list:
            d[v] = k
    return d


# From https://stackoverflow.com/a/26496899
def defaultdict_to_dict(d):
    if isinstance(d, defaultdict):
        d = {k: defaultdict_to_dict(v) for k, v in d.items()}
    return d


def sorted_nicely(l):
    """ Sort the given iterable in the way that humans expect.
        From http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def join_english_comma(array, last_join_word='and'):
    string = ''
    if array:
        if len(array) > 1:
            string += ', '.join(array[:-1])
            string += f" {last_join_word} {array[-1]}"
        else:
            string = array[0]
    return string


def full_class_name(klass):
    return klass.__module__ + '.' + klass.__name__


def import_class(full_class_path):
    modulename, classname = full_class_path.rsplit('.', 1)
    #logging.debug("import_class: modulename = %s", modulename)
    mod = importlib.import_module(modulename)
    return getattr(mod, classname)


def local_date_string() -> str:
    """ Returns eg '2022-07-18' """
    return localtime(timezone.now()).strftime("%Y-%m-%d")


class Struct:

    def __init__(self, **entries):
        self.__dict__.update(entries)


def time_since(start):
    end = time.time()
    return end - start


def get_and_log_time_since(start, name=''):
    ts = time_since(start)
    if name:
        name += ': '
    logging.info("%sTotal time taken: %.2f seconds", name, ts)
    return ts


def upper(string):
    if string:
        string = str(string).upper()
    return string


def is_url(url):
    parse_result = urlparse(url)
    return bool(parse_result.scheme in ('http', 'https') and parse_result.netloc)


def rgb_hex_to_tuples(rgb):
    rgb = rgb.replace('#', '')
    return bytes.fromhex(rgb)


def rgb_to_hex(red, green, blue):
    return "#%02x%02x%02x" % (red, green, blue)


def rgb_invert(rgb):
    red, green, blue = rgb_hex_to_tuples(rgb)
    inverted = (255 - red, 255 - green, 255 - blue)
    return rgb_to_hex(*inverted)


def get_all_subclasses(cls):
    """ From https://stackoverflow.com/a/17246726 """
    all_subclasses = set()

    for subclass in cls.__subclasses__():
        all_subclasses.add(subclass)
        all_subclasses.update(get_all_subclasses(subclass))

    return all_subclasses


def datetime_string_to_date(s):
    return parser.parse(s).date()


def none_to_blank_string(s):
    return s or ''


def calculate_age(born, died=None):
    """ https://stackoverflow.com/a/9754466/295724 """
    age = None
    if born:
        if died is None:
            today = date.today()
            age = today.year - born.year - ((today.month, today.day) < (born.month, born.day))
        else:
            age = died.year - born.year - ((died.month, died.day) < (born.month, born.day))
    return age


def all_equal(iterable):
    it = iter(iterable)
    start = next(it, None)
    return all(start == x for x in it)


def empty_dict():
    return dict()


T = TypeVar('T')
def empty_to_none(value: T) -> Optional[T]:
    """
    :param value: Any kind of value
    :return: None if value was whitespace, None, empty Dict or List, trimmed if value was a string, otherwise original value
    """

    # we wish to maintain False and 0 as a value
    if value is False:
        return False
    if value == 0:
        return 0

    if isinstance(value, str):
        value = value.strip()

    if not value:
        value = None

    return value


def pretty_label(label: str):
    label = label.replace('_', ' ')
    tidied = ''
    last_space = True
    for char in label:
        if last_space:
            char = char.upper()
            last_space = False
        if char == ' ':
            last_space = True
        tidied = tidied + char
    return tidied


def query_to_array(value):
    try:
        return json.loads(value)
    except:
        pass
    return [p.strip() for p in value.split(',')]


def is_not_none(obj):
    return obj is not None


def first(obj):
    if isinstance(obj, list):
        if len(obj) >= 1:
            return obj[0]
        return None
    return obj


def get_single_element(sequence):
    """ single element from set or list """
    size = len(sequence)
    if size != 1:
        msg = f"Expected single element sequence, but had {size} entries ({sequence})"
        raise ValueError(msg)
    return next(iter(sequence))


def nest_dict(flat_dict: dict) -> dict:
    nested = {}
    for full_path, value in flat_dict.items():
        path = full_path.split('.')
        nest_me = nested

        for idx, part in enumerate(path):
            if idx == len(path) - 1:
                nest_me[part] = value
            else:
                sub_dict = nest_me.get(part, {})
                nest_me[part] = sub_dict
                nest_me = sub_dict
    return nested


def iter_fixed_chunks(iterable, chunk_size):
    """ https://stackoverflow.com/a/22045226 """
    it = iter(iterable)
    return iter(lambda: tuple(islice(it, chunk_size)), ())


def nice_class_name(obj_or_klass):
    if isinstance(obj_or_klass, type):
        klass = obj_or_klass
    else:
        klass = obj_or_klass.__class__
    return klass.__name__


def single_quote(s):
    return f"'{s}'"


def double_quote(s):
    return f'"{s}"'


def filename_safe(filename) -> str:
    keepcharacters = {'.', '_'}
    filename = filename.replace(' ', '_')  # you can never trust spaces
    # leave room for an extension so make sure the filename is 250 characters
    filename = "".join(c for c in filename if c.isalnum() or c in keepcharacters).strip().lower()[0:250]
    return filename


def html_id_safe(text: str) -> str:
    """
    Makes a string that can be made into an HTML id and referenced easily
    """
    if not text:
        return str(uuid.uuid4())
    text = re.sub("[^0-9A-Za-z]", "-", text)
    text = re.sub("_{2,}", "-", text)
    if not text[0].isalpha():
        text = "x" + text
    return text


def html_link(url: str, title: str) -> SafeString:
    if not url:
        return mark_safe(title)
    return mark_safe(f"<a href='{url}'>{html.escape(title)}</a>")


T = TypeVar("T")


def batch_iterator(iterable: Iterable[T], batch_size: int = 10) -> Iterator[List[T]]:
    batch: List[T] = list()
    for record in iterable:
        batch.append(record)
        if len(batch) >= batch_size:
            yield batch
            batch = list()
    if batch:
        yield batch


def make_json_safe_in_place(obj):
    """
    converts all Decimals in a nested dict/array to floats.
    Specifically because ijson may return some
    """
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(value, Decimal):
                obj[key] = float(value)
            else:
                make_json_safe_in_place(value)
    elif isinstance(obj, list):
        for index, value in enumerate(obj):
            if isinstance(value, Decimal):
                obj[index] = float(value)
            else:
                make_json_safe_in_place(value)
    else:
        pass


class ChoicesEnum(Enum):
    """
    Extend this for enums that are going to work as Choices for models
    """

    @classmethod
    def choices(cls):
        return [(tag.value, tag.name) for tag in cls]


class ModelUtilsMixin:
    """
    Allows you to deal with an instance of a model or the key of the model
    """

    @classmethod
    def get(cls, value):
        if value is None:
            return None
        if isinstance(value, cls):
            return value
        if isinstance(value, (str, int)):
            return cls.objects.get(pk=value)
        raise ValueError(f'Expected {cls.__name__} or str or int, got {value}')


def execute_cmd(cmd: list, **kwargs) -> Tuple[int, Optional[str], Optional[str]]:
    shell = kwargs.get("shell", settings.POPEN_SHELL)

    if shell:
        command = ' '.join(cmd)
        logging.info('About to call %s', command)
        pipes = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
        logging.info('Completed')
    else:
        pipes = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)

    std_out, std_err = pipes.communicate()
    return pipes.returncode, std_out.decode() if std_out else None, std_err.decode() if std_err else None


Key = TypeVar("Key")
Value = TypeVar("Value")
Data = TypeVar("Data")

@dataclass
class Group(Generic[Key, Value]):
    key: Key
    values: Set[Value]

    @property
    def values_list(self):
        return sorted(self.values)

    def __lt__(self, other):
        return self.key < other.key


def group_data(data: Iterable[Data], key_func: Callable[[Data], Tuple[Key, Value]]) -> List[Group[Key, Value]]:
    group_dict = defaultdict(set)
    for element in data:
        key, value = key_func(element)
        group_dict[key].add(value)
    flat: List[Group[Key, Value]] = list()
    for key, values in group_dict.items():
        flat.append(Group(key, values))
    return flat


def group_by_key(qs: Iterable, key: attrgetter) -> Iterator[Tuple[Any, List]]:
    """
    Provide a sorted iterable value (like a queryset) and an attrgetter for getting values from it.
    Each time the attrgetter returns a value different from before, it will return a list
    of values where attrgetter returned the same.
    Note that if qs is not sorted you could get the same key value multiple times (just not right after each other)
    :param qs:
    :param key:
    :return:
    """
    last_value = None
    records = list()
    for record in qs:
        value = key(record)
        if value is None:
            raise ValueError('Grouping key must always have a value')
        if last_value is None or last_value != value:
            if records:
                yield last_value, records
            last_value = value
            records = list()
        records.append(record)
    if records:
        yield last_value, records


# note this tags expected in a single line of text
# don't catch too many tags in case you get some false positives
EXPECTED_HTML_TAGS_SINGLE_LINE = {'div', 'b', 'i', 'u', 'strong', 'em', 'sup', 'sub'}


def cautious_attempt_html_to_text(text: str, whitelist: Set[str] = None) -> str:
    """
    Given some text, and an expected whitelist of tags, will convert the text from possible HTML content to plain text.
    Converts things like &#039; to ' and will strip out expected tags, if a tag is found that's not in the whitelist
    the text will be returned untouched. This is to avoid treating text like "Patient had ouchies <painful>" as HTML
    :param text: Text that may contain HTML elements
    :param whitelist: HTML tags that we expect and want stripped out
    :return: Text stripped of HTML (if only tags present were whitelisted)
    """
    if whitelist is None:
        whitelist = EXPECTED_HTML_TAGS_SINGLE_LINE

    if not text:
        return text
    bs = BeautifulSoup(text, features="html.parser")
    for tag in bs.find_all():
        if tag.name.lower() not in whitelist:
            return text
    return bs.get_text()


def html_to_text(html: str, preserve_lines: bool = False) -> Optional[str]:
    if not html:
        return None
    bs = BeautifulSoup(f'<body>{html}</body>', features="html.parser")

    if not preserve_lines:
        return bs.get_text()
    else:

        def replace_with_newlines(element):
            text = ''
            for elem in element.recursiveChildGenerator():
                if isinstance(elem, str):
                    text += elem.strip() + " "
                elif elem.name == 'br':
                    text += '\n'
            return text.strip()

        def get_plain_text(soup):
            plain_text = ''
            lines = soup.find("body")
            return replace_with_newlines(lines)

        return get_plain_text(bs).strip()


def timestamp_as_number_formatter(row: Dict, field: str):
    """
    Useful for JQGrids where we want to send down time as a timestamp so we can put it in the user's timezone
    :param row: JQGrid row object
    :param field: field to access row (row[field] should be a timestamp)
    :return: unix time
    """
    val = row[field]
    return val.timestamp() if val else None


def get_subclasses(cls):
    """returns all subclasses of argument, cls"""
    if issubclass(cls, type):
        subclasses = cls.__subclasses__(cls)
    else:
        subclasses = cls.__subclasses__()
    for subclass in subclasses:
        subclasses.extend(get_subclasses(subclass))
    return subclasses


def count(obj: Any) -> int:
    if obj is None:
        return 0
    if isinstance(obj, int):
        return obj
    if isinstance(obj, QuerySet):
        return obj.count()
    return len(obj)


trailing_zeros_strip = re.compile("(.*?[.][0-9]*?)(0+)$")


def format_percent(number, is_unit=False) -> str:
    if is_unit:
        number *= 100
    return f"{format_significant_digits(number)}%"


def format_significant_digits(a_number, sig_digits=3) -> str:
    if a_number == 0:
        return "0"
    rounded_number = round(a_number, sig_digits - int(math.floor(math.log10(abs(a_number)))) - 1)
    rounded_number_str = "{:.12f}".format(rounded_number)
    if match := trailing_zeros_strip.match(rounded_number_str):
        rounded_number_str = match.group(1)
        if rounded_number_str[-1] == '.':
            rounded_number_str = rounded_number_str[:-1]

    return rounded_number_str


def delimited_row(data: list, delimiter: str = ',') -> str:
    out = io.StringIO()
    writer = csv.writer(out, delimiter=delimiter)
    writer.writerow(data)
    return out.getvalue()


def clean_string(input_string: str) -> str:
    """ Removes non-printable characters, strips whitespace """
    return re.sub(f'[^{re.escape(string.printable)}]', '', input_string.strip())


T = TypeVar("T")


class IterableTransformer(Generic[T], Iterable[T]):
    """
    Given an iterable and a transformer, makes a new iterable that will lazily transform the elements
    """

    class _IteratorTransformer:
        def __init__(self, iterator: Iterator, transform):
            self.iterator = iterator
            self.transform = transform

        def __next__(self):
            return self.transform(next(self.iterator))

    def __init__(self, iterable: Iterable, transform):
        self.iterable = iterable
        self.transform = transform

    def __iter__(self) -> Iterator[T]:
        return self._IteratorTransformer(iter(self.iterable), self.transform)


class IteratableStitcher(Generic[T], Iterable[T]):
    """
    Given a list of SORTED iterables, will give you the smallest element from any of them
    """

    class _IteratorStitcher:

        class _CachedIterator:

            def __init__(self, iteratable):
                self.iterator = iter(iteratable)
                self.finished = False
                self.cache = None
                self.fetch_next()

            def fetch_next(self):
                if not self.finished:
                    try:
                        self.cache = next(self.iterator)
                    except StopIteration:
                        self.cache = None
                        self.finished = True

            @property
            def preview(self):
                return self.cache

        def __init__(self, iteratables: List[Iterable[T]], comparison):
            self.iterators = [self._CachedIterator(iterable) for iterable in iteratables]
            self.comparison = comparison

        def __next__(self):
            min_ci: Optional[IteratableStitcher._IteratorStitcher._CachedIterator] = None
            min_value: Any = None

            for ci in self.iterators:
                value = ci.preview
                if value is not None:
                    if min_value is None or self.comparison(value, min_value):
                        min_value = value
                        min_ci = ci

            if not min_ci:
                raise StopIteration()
            min_ci.fetch_next()
            return min_value

    def __init__(self, iterables: List[Iterable[T]], comparison=operator.__lt__):
        self.iterables = iterables
        self.comparison = comparison

    def __iter__(self) -> Iterator[T]:
        return self._IteratorStitcher(iteratables=self.iterables, comparison=self.comparison)


class Constant:
    """ Used for creating non-enum properties in Enums
        From https://stackoverflow.com/a/18035135 """
    def __init__(self, value):
        self.value = value
    def __get__(self, *args):
        return self.value
    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.value)


class ArrayLength(models.Func):
    """
    Can annotate with array length now e.g.
    MyModel.objects.all().annotate(field_len=ArrayLength('field')).order_by('field_len')
    """
    function = 'CARDINALITY'


@dataclass
class DebugTime:
    description: str
    durations: List[timedelta] = field(default_factory=list)
    occurrences: int = 0

    def tick(self, duration: timedelta):
        self.durations.append(duration)
        self.occurrences += 1

    @property
    def duration(self):
        return reduce(lambda x, total: x + total, self.durations, timedelta())

    def __str__(self):
        if self.occurrences == 1:
            return f"{self.duration} - {self.description}"
        else:
            return f"{self.duration / self.occurrences} (x {self.occurrences}) - {self.description} "


class DebugTimer:

    def __init__(self):
        self.start = datetime.now()
        self.times: Dict[str, DebugTime] = dict()

    def tick(self, description: str):
        now = datetime.now()
        duration = now - self.start

        debug_time: DebugTime
        if existing := self.times.get(description):
            debug_time = existing
        else:
            debug_time = DebugTime(description)
            self.times[description] = debug_time

        debug_time.tick(duration)
        self.start = now

    def __str__(self):
        return "\n".join((str(debug_time) for debug_time in self.times.values()))


class NullTimer(DebugTimer):

    def tick(self, description: str):
        pass

    def __str__(self):
        return "NullTimer"


DebugTimer.NullTimer = NullTimer()


class LimitedCollection:

    def __init__(self, data: List[Any], limit: int):
        if data is None:
            data = []
        self.true_count = len(data)
        self.data = data
        if len(data) > limit:
            self.data = data[0:limit]
        self.limit = limit

    @property
    def limit_str(self):
        return f"Limiting results to {len(self.data)} of {self.true_count}"

    @property
    def is_limited(self) -> bool:
        return len(self.data) != self.true_count

    def __len__(self):
        return self.true_count

    def __iter__(self):
        return iter(self.data)

    def __getitem__(self, item):
        return self.data[item]

    def __bool__(self):
        return self.true_count > 0


class LazyAttribute:
    """
    LazyAttributes allow you to easily make a context map for a view where properties
    will only be evaluated when accessed.
    Importantly this works for attributes marked with @lazy which are ironically eagerly
    evaluated by SimpleLazyObject

    class SomeClass:
        @lazy
        def lazy_prop():
            return "foo"

        def method():
            return "bar"

    instance = SomeClass()
    LazySimpleObject(instance.lazy_prop) # is eager
    LazySimpleObject(lambda: instance.lazy_prop) # is eager
    LazySimpleObject(instance.method) # is lazy
    LazySimpleObject(LazyAttribute(instance, "lazy_prop") # is lazy
    """

    def __init__(self, obj: Any, attr: str):
        self.obj = obj
        self.attr = attr

    def eval(self):
        attr = getattr(self.obj, self.attr)
        if callable(attr):
            attr = attr()
        return attr

    @staticmethod
    def lazy_context(obj: Any, attributes: List[str]) -> Dict[str, SimpleLazyObject]:
        """
        Create a dict used for contexts retrieve attributes from obj and having them only
        evaluate when the template first mentions them
        """
        all_attributes = set(dir(obj))

        context = {}
        for attribute in attributes:
            if attribute not in all_attributes:
                raise ValueError(f"{obj} does not have attribute {attribute}")
            lazy_att = LazyAttribute(obj, attribute)
            context[attribute] = SimpleLazyObject(lazy_att.eval)
        return context


P = TypeVar('P')


def segment(iterable: Iterable[P], filter: Callable[[P], bool]) -> Tuple[List[P], List[P]]:
    """
    :param iterable An iterable bunch of data to be split into two
    :param filter A filter to run over each element of iterable, to put it into a pass or fail list
    :returns two lists, the first being elements that passed the filter, the second being ones that failed
    """
    passes: List[P] = list()
    fails: List[P] = list()
    for element in iterable:
        if filter(element):
            passes.append(element)
        else:
            fails.append(element)
    return passes, fails


def export_column(label: Optional[str] = None, sub_data: Optional[Type] = None, categories: Dict[str, Any] = None):
    """
    Extend ExportRow and annotate methods with export_column.
    The order of defined methods determines the order that the results will appear in an export file
    :param label: The label that will appear in the CSV header (defaults to method name if not provided)
    :param sub_data: An optional SubType of another ExportRow for nested data
    """

    def decorator(method):
        def wrapper(*args, **kwargs):
            return method(*args, **kwargs)
        # have to cache the line number of the source method, otherwise we just get the line number of this wrapper
        wrapper.line_number = inspect.getsourcelines(method)[1]
        wrapper.label = label or method.__name__
        wrapper.__name__ = method.__name__
        wrapper.is_export = True
        wrapper.categories = categories
        wrapper.sub_data = sub_data
        return wrapper
    return decorator


class ExportRow:

    @staticmethod
    def get_export_methods(klass, categories: Optional[Dict[str, Any]] = None):
        if not hasattr(klass, 'export_methods'):
            export_methods = [func for _, func in inspect.getmembers(klass, lambda x: getattr(x, 'is_export', False))]
            export_methods.sort(key=lambda x: x.line_number)
            klass.export_methods = export_methods

            if not klass.export_methods:
                raise ValueError(f"ExportRow class {klass} has no @export_columns, make sure you annotate @export_column(), not @export_column")

        export_methods = klass.export_methods
        if categories:
            def passes_filter(export_method) -> bool:
                nonlocal categories
                export_categories = export_method.categories or dict()
                for key, value in categories.items():
                    if not export_categories.get(key) == value:
                        return False
                return True

            export_methods = [em for em in export_methods if passes_filter(em)]

        return export_methods

    @classmethod
    def _data_generator(cls, data: Iterable[Any]) -> Iterator[Any]:
        for row_data in data:
            if row_data is None:
                continue
            if not isinstance(row_data, cls):
                # it's expected that the class can be initiated with each "row" in data
                row_data = cls(row_data)
            yield row_data

    @classmethod
    def csv_generator(cls, data: Iterable[Any], delimiter=',', include_header=True, categories: Optional[Dict[str, Any]] = None) -> Iterator[str]:
        try:
            if include_header:
                yield delimited_row(cls.csv_header(categories=categories), delimiter=delimiter)
            for row_data in cls._data_generator(data):
                yield delimited_row(row_data.to_csv(categories=categories), delimiter=delimiter)
        except:
            from library.log_utils import report_exc_info
            report_exc_info(extra_data={"activity": "Exporting"})
            yield "** File terminated due to error"
            raise

    @classmethod
    def json_generator(cls, data: Iterable[Any], records_key: str = "records") -> Iterator[str]:
        first_row = True
        try:
            yield f'{{"{records_key}": ['
            for row_data in cls._data_generator(data):
                yield (', ' if not first_row else '') + json.dumps(row_data.to_json())
                first_row = False
            yield ']}}'
        except:
            from library.log_utils import report_exc_info
            report_exc_info(extra_data={"activity": "Exporting"})
            yield "\"error\"** File terminated due to error"
            raise

    @classmethod
    def csv_header(cls, categories: Optional[Dict[str, Any]] = None) -> List[str]:
        row = list()
        for method in ExportRow.get_export_methods(cls, categories=categories):
            label = method.label or method.__name__
            if sub_data := method.sub_data:
                sub_header = sub_data.csv_header(categories=categories)
                for sub in sub_header:
                    row.append(label + "." + sub)
            else:
                row.append(label)
        return row

    def to_csv(self, categories: Optional[Dict[str, Any]] = None) -> List[str]:
        row = list()
        for method in ExportRow.get_export_methods(self.__class__, categories=categories):
            result = method(self)
            if sub_data := method.sub_data:
                if result is None:
                    for entry in sub_data.csv_header(categories=categories):
                        row.append("")
                else:
                    row += result.to_csv()
            else:
                row.append(result)
        return row

    def to_json(self, categories: Optional[Dict[str, Any]] = None) -> JsonObjType:
        row = dict()
        for method in ExportRow.get_export_methods(self.__class__, categories=categories):
            result = method(self)
            value: Any
            if result is None:
                value = None
            elif method.sub_data:
                value = result.to_json(categories=categories)
            else:
                value = result

            if value == "":
                value = None
            row[method.__name__] = value

        return row

    @classmethod
    def streaming(cls, request: HttpRequest, data: Iterable[Any], filename: str, categories: Optional[Dict[str, Any]] = None):
        if request.GET.get('format') == 'json':
            return cls.streaming_json(data, filename, categories=categories)
        else:
            return cls.streaming_csv(data, filename, categories=categories)

    @classmethod
    def streaming_csv(cls, data: Iterable[Any], filename: str, categories: Optional[Dict[str, Any]] = None):
        date_str = local_date_string()

        response = StreamingHttpResponse(cls.csv_generator(data, categories=categories), content_type='text/csv')
        response['Content-Disposition'] = f'attachment; filename="{filename}_{settings.SITE_NAME}_{date_str}.csv"'
        return response

    @classmethod
    def streaming_json(cls, data: Iterable[Any], filename: str, records_key: str = None, categories: Optional[Dict[str, Any]] = None):
        date_str = local_date_string()

        if not records_key:
            records_key = filename.replace(" ", "_")

        response = StreamingHttpResponse(cls.json_generator(data, records_key, categories=categories), content_type='application/json')
        response['Content-Disposition'] = f'attachment; filename="{filename}_{settings.SITE_NAME}_{date_str}.json"'
        return response


class VarsDict:

    def get(self, key, default=None):
        return vars(self).get(key, default)

    def __contains__(self, item):
        return bool(self.get(item))

    def __getitem__(self, item):
        return vars(self)[item]


@dataclass(frozen=True)
class DiffTextSegment:
    operation: str
    text: str

    @property
    def operation_name(self):
        if self.operation == ' ':
            return 'same'
        elif self.operation == '-':
            return 'subtract'
        elif self.operation == '+':
            return 'add'
        else:
            return self.operation

    def __str__(self):
        return f"{self.operation} {self.text}"


class DiffBuilder:

    def __init__(self):
        self.diff_segments: List[DiffTextSegment] = list()
        self.same_text = ''
        self.add_text = ''
        self.sub_text = ''

    def apply(self):
        if sub_text := self.sub_text:
            self.diff_segments.append(DiffTextSegment(operation='-', text=sub_text))
        if add_text := self.add_text:
            self.diff_segments.append(DiffTextSegment(operation='+', text=add_text))
        if same_text := self.same_text:
            self.diff_segments.append(DiffTextSegment(operation=' ', text=same_text))
        self.sub_text = ''
        self.add_text = ''
        self.same_text = ''

    def _has_pending_add_subtract(self):
        return bool(self.add_text) or bool(self.sub_text)

    def _has_pending_same(self):
        return bool(self.same_text)

    def append(self, diff_text: str):
        operation = diff_text[0]
        text = diff_text[2:]
        if operation == '?':
            # '? 'line not present in either input sequence
            # Lines beginning with ‘?’ attempt to guide the eye to intraline differences,
            # and were not present in either input sequence. These lines can be confusing if the sequences contain tab characters.
            return
        if operation == ' ':
            if self._has_pending_add_subtract():
                if not text.strip():
                    # if we have "bear down" -> "fat dog" it would result in -bear +fat (same space) -down +dog
                    # much cleaner to have that as -"bear down" +"fat dog"
                    # might need to look into having whitespace_text that we track separately
                    self.sub_text += text
                    self.add_text += text
                    return

                self.apply()
            self.same_text += text
        else:
            if self._has_pending_same():
                self.apply()
            if operation == '-':
                self.sub_text += text
            else:
                self.add_text += text

    def __iter__(self):
        return iter(self.diff_segments)

    def __bool__(self):
        return bool(self.diff_segments)


def diff_text(a: str, b: str) -> DiffBuilder:

    def _tokenize(text: str) -> List[str]:
        return re.split(r'(\s)', text)

    diff_builder = DiffBuilder()
    for diff_chars in difflib.Differ().compare(_tokenize(a), _tokenize(b)):
        diff_builder.append(diff_chars)
    diff_builder.apply()
    print(diff_builder.diff_segments)
    return diff_builder
