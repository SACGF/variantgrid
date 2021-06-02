import csv
import io
import operator
import math
from collections import defaultdict
from datetime import date, datetime
from operator import attrgetter
from urllib.parse import urlparse

from bs4 import BeautifulSoup
from dateutil import parser
from decimal import Decimal
from django.core.serializers import serialize
from django.db.models.query import QuerySet
from django.db import models
from enum import Enum
from itertools import islice
from json.encoder import JSONEncoder
from typing import TypeVar, Optional, Iterator, Tuple, Any, List, Iterable, Set, Dict, Union
import hashlib
import importlib
import json
import logging
import re
import subprocess
import time
from django.conf import settings
from django.utils import html
from django.utils.functional import SimpleLazyObject
from django.utils.safestring import SafeString, mark_safe

FLOAT_REGEX = r'([-+]?[0-9]*\.?[0-9]+.|Infinity)'


# isn't restrictive enough, but helps documention
JSON = Union[Dict[str,Any], List[Any], int, float, str]


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


def html_link(url: str, title: str) -> SafeString:
    if not url:
        return mark_safe(title)
    return mark_safe(f"<a href='{url}'>{html.escape(title)}</a>")


def batch_iterator(iterable, batch_size: int = 10):
    batch = []
    for record in iterable:
        batch.append(record)
        if len(batch) >= batch_size:
            yield batch
            batch = []
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
EXPECTED_HTML_TAGS_SINGLE_LINE = set([
    'div', 'b', 'i', 'u', 'strong', 'em'
])


def cautious_attempt_html_to_text(text: str, whitelist: Set[str] = None) -> str:
    """
    Given some text, and an expected whitelist of tags, will convert the text from possible HTML content to plain text.
    Converts things like &#039; to ' and will strip out execpted tags, if a tag is found that's not in the whitelist
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


class IterableTransformer:
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

    def __iter__(self):
        return self._IteratorTransformer(iter(self.iterable), self.transform)


class IteratableStitcher:
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

        def __init__(self, iteratables: List[Iterable], comparison):
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

    def __init__(self, iterables: List[Iterable], comparison=operator.__lt__):
        self.iterables = iterables
        self.comparison = comparison

    def __iter__(self):
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


class DebugTimer:

    def __init__(self):
        self.start = datetime.now()

    def tick(self, description: str):
        now = datetime.now()
        duration = now - self.start
        print(f"{description} {duration}")
        self.start = now


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
        return self.limit

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
