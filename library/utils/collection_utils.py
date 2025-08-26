import abc
import itertools
import operator
import re
from abc import ABC
from collections import defaultdict
from dataclasses import dataclass, field
from html import escape
from itertools import islice
from typing import Iterable, Iterator, TypeVar, Any, Generic, Callable, Optional, Sequence, Union

from django.utils.functional import SimpleLazyObject
from django.utils.safestring import SafeString

DictKey = TypeVar("DictKey")
DictVal = TypeVar("DictVal")


def invert_dict(source_dict: dict[DictKey, DictVal]) -> dict[DictVal, DictKey]:
    return {v: k for k, v in source_dict.items()}


def invert_dict_of_lists(source_dict):
    d = {}
    for k, v_list in source_dict.items():
        for v in v_list:
            d[v] = k
    return d


# From https://stackoverflow.com/a/26496899
def defaultdict_to_dict(d: defaultdict) -> dict:
    if isinstance(d, defaultdict):
        d = {k: defaultdict_to_dict(v) for k, v in d.items()}
    return d


def sorted_nicely(l: Iterable) -> Sequence:
    """ Sort the given iterable in the way that humans expect.
        From http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html """
    def convert(text: str) -> Union[int, str]:
        return int(text) if text.isdigit() else text

    def alphanum_key(key: str) -> list[str]:
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(l, key=alphanum_key)


T = TypeVar("T")


def first(obj: Iterable[T]) -> T:
    # FIXME, use collection_tools.first()
    for element in obj:
        return element
    return None


def get_single_element(sequence: Sequence[T]) -> T:
    """ single element from set or list """
    # FIXME this is probably better done by just trying next( ) and catching and re-throwing the exception if sie is 0
    size = len(sequence)
    if size != 1:
        msg = f"Expected single element sequence, but had {size} entries ({sequence})"
        raise ValueError(msg)
    return next(iter(sequence))


def nest_dict(flat_dict: dict[DictKey, DictVal]) -> dict[DictKey, DictVal]:
    """
    :param flat_dict: A dictionary where all the keys are in the format of "a.b": x, "a.c": y
    :return: A nested dictionary e.g. {"a": {"b": x}, {"c": y}}
    """
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


def iter_fixed_chunks(iterable: Iterable[T], chunk_size: int) -> Iterator[tuple[T]]:
    """ https://stackoverflow.com/a/22045226 """
    it = iter(iterable)
    return iter(lambda: tuple(islice(it, chunk_size)), ())


def batch_iterator(iterable: Iterable[T], batch_size: int = 10) -> Iterator[list[T]]:
    batch: list[T] = []
    for record in iterable:
        batch.append(record)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


Key = TypeVar("Key")
Value = TypeVar("Value")
Data = TypeVar("Data")


@dataclass
class Group(Generic[Key, Value]):
    key: Key
    values: set[Value]

    @property
    def values_list(self):
        return sorted(self.values)

    def __lt__(self, other):
        return self.key < other.key


def group_data(data: Iterable[Data], key_func: Callable[[Data], tuple[Key, Value]]) -> list[Group[Key, Value]]:
    # deprecated, use itertools groupby
    group_dict = defaultdict(set)
    for element in data:
        key, value = key_func(element)
        group_dict[key].add(value)
    flat: list[Group[Key, Value]] = []
    for key, values in group_dict.items():
        flat.append(Group(key, values))
    return flat


def group_by_key(qs: Iterable, key: Callable) -> Iterator[tuple[Any, list]]:
    # DEPRECATED - use itertoool.group_by instead
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
    records = []
    for record in qs:
        value = key(record)
        if value is None:
            raise ValueError('Grouping key must always have a value')
        if last_value is None or last_value != value:
            if records:
                yield last_value, records
            last_value = value
            records = []
        records.append(record)
    if records:
        yield last_value, records


def sort_and_group_by(collection, key_func: Callable[[Any], Any]) -> Iterable[tuple[Any, list]]:
    sorted_data = sorted(collection, key=key_func)
    return itertools.groupby(sorted_data, key=key_func)


def all_equal(iterable):
    it = iter(iterable)
    start = next(it, None)
    return all(start == x for x in it)


TransformInput = TypeVar("TransformInput")


class IterableTransformer(Generic[T], Iterable[T]):
    """
    Given an iterable and a transformer, makes a new iterable that will lazily transform the elements
    """

    class _IteratorTransformer(Iterator[T]):
        def __init__(self, iterator: Iterator, transform):
            self.iterator = iterator
            self.transform = transform

        def __next__(self):
            return self.transform(next(self.iterator))

    def __init__(self, iterable: Iterable[TransformInput], transform: Callable[[TransformInput], T]):
        self.iterable = iterable
        self.transform = transform

    def __iter__(self) -> Iterator[T]:
        return self._IteratorTransformer(iter(self.iterable), self.transform)


class IterableStitcher(Generic[T], Iterable[T]):
    """
    Given a list of SORTED iterables, will give you the smallest element from any of them
    """

    class _IteratorStitcher(Iterator[T]):

        class _CachedIterator:

            def __init__(self, iterable):
                self.iterator = iter(iterable)
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

        def __init__(self, iterables: list[Iterable[T]], comparison: Callable[[T, T], bool]):
            self.iterators = [self._CachedIterator(iterable) for iterable in iterables]
            self.comparison = comparison

        def __next__(self):
            min_ci: Optional[IterableStitcher._IteratorStitcher._CachedIterator] = None
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

    def __init__(self, iterables: list[Iterable[T]], comparison=operator.__lt__):
        self.iterables = iterables
        self.comparison = comparison

    def __iter__(self) -> Iterator[T]:
        return self._IteratorStitcher(iterables=self.iterables, comparison=self.comparison)


class LimitedCollection(Generic[T]):

    def __init__(self, data: list[T], limit: int):
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

    def __iter__(self) -> Iterator[T]:
        return iter(self.data)

    def __getitem__(self, item) -> T:
        return self.data[item]

    def __bool__(self):
        return self.true_count > 0


def segment(iterable: Iterable[T], filter_func: Callable[[T], bool]) -> tuple[list[T], list[T]]:
    """
    :param iterable An iterable bunch of data to be split into two
    :param filter_func A filter to run over each element of iterable, to put it into a pass or fail list
    :returns two lists, the first being elements that passed the filter, the second being ones that failed
    """
    passes: list[T] = []
    fails: list[T] = []
    for element in iterable:
        if filter_func(element):
            passes.append(element)
        else:
            fails.append(element)
    return passes, fails


def flatten_nested_lists(iterable) -> list:
    # collapses lists of lists, and filters out Nones
    def _flatten_generator(flatten_me):
        for item in flatten_me:
            if isinstance(item, list):
                yield from _flatten_generator(item)
            elif item is not None:
                yield item

    return list(_flatten_generator(iterable))


ListItem = TypeVar("ListItem")


def remove_duplicates_from_list(source_list: list[ListItem]) -> list[ListItem]:
    seen: set[ListItem] = set()
    output: list[ListItem] = []
    for item in source_list:
        if item not in seen:
            seen.add(item)
            output.append(item)
    return output


class VarsDict:

    def get(self, key, default=None):
        return vars(self).get(key, default)

    def __contains__(self, item):
        return bool(self.get(item))

    def __getitem__(self, item):
        return vars(self)[item]


class LazyAttribute:
    """
    LazyAttributes allow you to easily make a context map for a view where properties
    will only be evaluated when accessed.
    Importantly this works for attributes marked with @cached_property which are ironically eagerly
    evaluated by SimpleLazyObject

    class SomeClass:
        @cached_property
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
    def lazy_context(obj: Any, attributes: list[str]) -> dict[str, SimpleLazyObject]:
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


class FormerTuple(ABC):
    """
    Base class for converting tuples (or named tuples) to data classes while maintaining decomposition
    e.g.
    @dataclass(frozen=True)
    class VariantCoordinate(FormerTuple):
        chrom: str
        start: int
        end: int
        alt: str
        ref: str

        @property
        def as_tuple():
            return (self.chrom, self.start, self.end, self.alt, self.ref)

    # existing code (still works)
    chrom, end, alt, ref, svlen = some_variant_coordinate
    """

    @property
    @abc.abstractmethod
    def as_tuple(self) -> tuple:
        pass

    def __iter__(self):
        return iter(self.as_tuple)

    def __getitem__(self, item):
        return self.as_tuple[item]

    def __eq__(self, other):
        # default implementation assumes all parts implement __eq__
        if not type(self) is type(other):
            return False
        return self.as_tuple == other.as_tuple


@dataclass(frozen=True)
class RowSpanCellValue:
    value: Any = None
    css_classes: list[str] = field(default_factory=list)
    href: Optional[str] = None


@dataclass
class RowSpanCell:
    repeat_of: Optional['RowSpanCell'] = None
    cell: RowSpanCellValue = None
    row_span: int = 1
    start_tr: bool = False
    end_tr: bool = False

    def add_row_span(self):
        if repeat_of := self.repeat_of:
            repeat_of.add_row_span()
        else:
            self.row_span += 1

    @property
    def link_contents(self):
        if self.href:
            href_escaped = self.href.replace("\"", "'")
            if self.href.startswith("http") or self.href.startswith("/"):
                return SafeString(f""" class="modal-link" data-toggle="ajax-modal" data-size="lg" data-title="Comments" data-href="{href_escaped}" """)
            else:
                return SafeString(f"onclick=\"{href_escaped}\" class=\"filter-link\" ")

    @property
    def value(self):
        return self.cell.value

    @property
    def css_class(self):
        if self.cell.css_classes:
            return " ".join(self.cell.css_classes)
        return ""

    @property
    def href(self):
        return self.cell.href


class RowSpanTable:

    def __init__(self, column_count: int):
        self.column_count = column_count
        self.rows = []
        self.last_row = None

    def add_cell(self, column: int, cell: RowSpanCellValue):
        end_tr = column == self.column_count - 1
        if self.last_row is None:
            self.last_row = [None] * self.column_count
            self.last_row[column] = RowSpanCell(cell=cell, start_tr=True, end_tr=end_tr)
            self.rows.append(self.last_row)
            return

        if self.last_row[column] is None:
            self.last_row[column] = RowSpanCell(cell=cell, end_tr=end_tr)
        else:
            new_row = [None] * self.column_count
            for i in range(0, column):
                target = RowSpanCell(repeat_of=self.last_row[i])
                target.add_row_span()
                new_row[i] = target

            new_row[column] = RowSpanCell(cell=cell, start_tr=True, end_tr=end_tr)
            self.last_row = new_row
            self.rows.append(new_row)

    def cells(self):
        for row in self.rows:
            for cell in row:
                if cell.repeat_of is None:
                    yield cell

