import functools
import gzip
import io
import itertools
import json
import time
from enum import Enum
from json.encoder import JSONEncoder
from typing import TypeVar, Optional, Any
from urllib.parse import urlparse

from django.core.serializers import serialize
from django.db.models.query import QuerySet
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

FLOAT_REGEX = r'([-+]?[0-9]*\.?[0-9]+.|Infinity)'


class DjangoJSONEncoder(JSONEncoder):

    def default(self, o):
        if isinstance(o, QuerySet):
            # `default` must return a python serializable
            # structure, the easiest way is to load the JSON
            # string produced by `serialize` and return it
            return json.loads(serialize('json', o))
        return JSONEncoder.default(self, o)


class Struct:

    def __init__(self, **entries):
        self.__dict__.update(entries)


def is_url(url: str) -> bool:
    parse_result = urlparse(url)
    return bool(parse_result.scheme in ('http', 'https') and parse_result.netloc)


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


def is_not_none(obj):
    return obj is not None


class ChoicesEnum(Enum):
    """
    Deprecated: Use models.TextChoices now instead
    Extend this for enums that are going to work as Choices for models
    """

    @classmethod
    def choices(cls):
        return [(tag.value, tag.name) for tag in cls]


def count(obj: Any) -> int:
    if obj is None:
        return 0
    if isinstance(obj, int):
        return obj
    if isinstance(obj, QuerySet):
        return obj.count()
    return len(obj)


class Constant:
    """ Used for creating non-enum properties in Enums
        From https://stackoverflow.com/a/18035135 """

    def __init__(self, value):
        self.value = value

    def __get__(self, *args):
        return self.value

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self.value)


class WrappablePartial(functools.partial):
    """
    functools.partial doesn't work great on decorated methods (complains about lack of __module__)
    this extension of partial fixes that.
    Partial is having a method(arg1, arg2) and turning it into method(arg2) - where arg1 is a constant
    """

    @property
    def __module__(self):
        return self.func.__module__

    @property
    def __name__(self):
        return "functools.partial({}, *{}, **{})".format(
            self.func.__name__,
            self.args,
            self.keywords
        )

    @property
    def __doc__(self):
        return self.func.__doc__


def invalidate_cached_property(obj, property_name: str):
    """ Invalidates without throwing exception if not defined """
    if property_name in obj.__dict__:
        delattr(obj, property_name)


def update_dict_of_dict_values(dict_to_update: dict[Any, dict], new_values: dict[Any, dict]):
    for k, v in new_values.items():
        old_values = dict_to_update.get(k, {})
        old_values.update(v)
        dict_to_update[k] = old_values


_retry = Retry(total=5, read=5, backoff_factor=0, allowed_methods=["GET"])
_ses = Session()
_ses.mount("http://", HTTPAdapter(max_retries=_retry))
_ses.mount("https://", HTTPAdapter(max_retries=_retry))

def iter_http_lines(url, timeout=60, attempts=5, backoff=2, encoding="utf-8", session=_ses):
    for n in itertools.count():
        try:
            with session.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                binary = r.raw
                if r.headers.get("Content-Encoding") == "gzip":
                    binary.decode_content = True
                elif url.endswith(".gz"):
                    binary = gzip.GzipFile(fileobj=binary)
                text = io.TextIOWrapper(binary, encoding=encoding)
                for line in text:
                    yield line.rstrip("\n")
            return
        except Exception:
            if n + 1 >= attempts:
                raise
            time.sleep(backoff * 2**n)
