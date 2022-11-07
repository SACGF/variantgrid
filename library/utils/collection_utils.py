import operator
import re
from collections import defaultdict
from dataclasses import dataclass
from itertools import islice
from typing import Iterable, Iterator, List, TypeVar, Any, Generic, Set, Callable, Tuple, Optional, Dict, Sequence, \
    Union

from django.utils.functional import SimpleLazyObject

DictKey = TypeVar("DictKey")
DictVal = TypeVar("DictVal")


def invert_dict(source_dict: Dict[DictKey, DictVal]) -> Dict[DictVal, DictKey]:
    return {v: k for k, v in source_dict.items()}


def invert_dict_of_lists(source_dict):
    d = {}
    for k, v_list in source_dict.items():
        for v in v_list:
            d[v] = k
    return d


# From https://stackoverflow.com/a/26496899
def defaultdict_to_dict(d: defaultdict) -> Dict:
    if isinstance(d, defaultdict):
        d = {k: defaultdict_to_dict(v) for k, v in d.items()}
    return d


def sorted_nicely(l: Iterable) -> Sequence:
    """ Sort the given iterable in the way that humans expect.
        From http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html """
    def convert(text: str) -> Union[int, str]:
        return int(text) if text.isdigit() else text

    def alphanum_key(key: str) -> List[str]:
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(l, key=alphanum_key)


def first(obj):
    # FIXME, seems inferior to get_single_element in every way
    # (note this will return None instead of ValueError)
    if isinstance(obj, list):
        if len(obj) >= 1:
            return obj[0]
        return None
    return obj


T = TypeVar("T")


def get_single_element(sequence: Sequence[T]) -> T:
    """ single element from set or list """
    # FIXME this is probably better done by just trying next( ) and catching and re-throwing the exception if sie is 0
    size = len(sequence)
    if size != 1:
        msg = f"Expected single element sequence, but had {size} entries ({sequence})"
        raise ValueError(msg)
    return next(iter(sequence))


def nest_dict(flat_dict: Dict[DictKey, DictVal]) -> Dict[DictKey, DictVal]:
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


def iter_fixed_chunks(iterable: Iterable[T], chunk_size: int) -> Iterator[List[T]]:
    """ https://stackoverflow.com/a/22045226 """
    it = iter(iterable)
    return iter(lambda: tuple(islice(it, chunk_size)), ())


def batch_iterator(iterable: Iterable[T], batch_size: int = 10) -> Iterator[List[T]]:
    batch: List[T] = []
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
    values: Set[Value]

    @property
    def values_list(self):
        return sorted(self.values)

    def __lt__(self, other):
        return self.key < other.key


def group_data(data: Iterable[Data], key_func: Callable[[Data], Tuple[Key, Value]]) -> List[Group[Key, Value]]:
    # deprecated, use itertools groupby
    group_dict = defaultdict(set)
    for element in data:
        key, value = key_func(element)
        group_dict[key].add(value)
    flat: List[Group[Key, Value]] = []
    for key, values in group_dict.items():
        flat.append(Group(key, values))
    return flat


def group_by_key(qs: Iterable, key: Callable) -> Iterator[Tuple[Any, List]]:
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


def all_equal(iterable):
    it = iter(iterable)
    start = next(it, None)
    return all(start == x for x in it)


def empty_dict() -> Dict:
    # If you want an empty_dict as a default function parameter
    return {}


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

        def __init__(self, iterables: List[Iterable[T]], comparison: Callable[[T, T], bool]):
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

    def __init__(self, iterables: List[Iterable[T]], comparison=operator.__lt__):
        self.iterables = iterables
        self.comparison = comparison

    def __iter__(self) -> Iterator[T]:
        return self._IteratorStitcher(iterables=self.iterables, comparison=self.comparison)


class LimitedCollection(Generic[T]):

    def __init__(self, data: List[T], limit: int):
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


def segment(iterable: Iterable[T], filter: Callable[[T], bool]) -> Tuple[List[T], List[T]]:
    """
    :param iterable An iterable bunch of data to be split into two
    :param filter A filter to run over each element of iterable, to put it into a pass or fail list
    :returns two lists, the first being elements that passed the filter, the second being ones that failed
    """
    passes: List[T] = []
    fails: List[T] = []
    for element in iterable:
        if filter(element):
            passes.append(element)
        else:
            fails.append(element)
    return passes, fails


def flatten_nested_lists(iterable) -> List:
    # collapses lists of lists, and filters out Nones
    def _flatten_generator(flatten_me):
        for item in flatten_me:
            if isinstance(item, list):
                for sub_item in _flatten_generator(item):
                    yield sub_item
            elif item is not None:
                yield item

    return list(_flatten_generator(iterable))


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
