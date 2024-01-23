from abc import abstractmethod, ABC
from dataclasses import dataclass
from typing import Any, Optional

from django.conf import settings
from django.db.models import Manager, QuerySet
from frozendict import frozendict
from threadlocals.threadlocals import set_request_variable, get_request_variable


@dataclass(frozen=True)
class CachedQuery:
    model: Any
    args: Any
    kwargs: Any

    def __str__(self):
        return f"{self.model} {self.args} {self.kwargs}"


_CACHED_QUERIES: dict[CachedQuery, Any] = {}
# TODO, could make this a LRU cache so it doesn't grow infinitely for larger models


class QuerySetCachingBase(QuerySet, ABC):

    @abstractmethod
    def store_in_cache(self, cq: CachedQuery, result: QuerySet):
        pass

    @abstractmethod
    def retrieve_from_cache(self, cq: CachedQuery):
        pass

    def get(self, *args, **kwargs):
        cq: Optional[CachedQuery] = None
        kwargs_fd = None
        if kwargs:
            # turn into a frozen dict so we can hash it
            kwargs_fd = frozendict(kwargs)
        try:
            cq = CachedQuery(model=self.model, args=args, kwargs=kwargs_fd)
        except TypeError:
            # un-hashable argument
            pass

        if cq:
            existing = self.retrieve_from_cache(cq)
            if existing is not None:
                return existing
            else:
                result = super().get(*args, **kwargs)
                self.store_in_cache(cq, result)
                return result
        else:
            return super().get(*args, **kwargs)


class QuerySetImmutable(QuerySetCachingBase):

    # def __init__(self, model=None, query=None, using=None, hints=None):
    #     super().__init__(model=model, query=query, using=using, hints=hints)

    def store_in_cache(self, cq: CachedQuery, result: QuerySet):
        _CACHED_QUERIES[cq] = result

    def retrieve_from_cache(self, cq: CachedQuery):
        return _CACHED_QUERIES.get(cq)


class QuerySetRequestCache(QuerySetCachingBase):

    def store_in_cache(self, cq: CachedQuery, result: QuerySet):
        try:
            set_request_variable(key=cq, val=result, use_threadlocal_if_no_request=False)
        except RuntimeError:
            pass

    def retrieve_from_cache(self, cq: CachedQuery):
        try:
            return get_request_variable(key=cq, default=None, use_threadlocal_if_no_request=False)
        except RuntimeError:
            return None


class ObjectManagerCachingImmutable(Manager):
    """
    Best used as both the default objects and as Meta._base_manager_name = 'objects'
    Will cache results for queries - so only do this on models that are read from but never written to.

    We could invalidate the cache when save() is called, but that only works if the save is done in the only instance
    that has already cached values.

    Note that currently it only caches the results of get(...)
    """

    def __init__(self):
        super().__init__()
        if not settings.UNIT_TEST:
            self._queryset_class = QuerySetImmutable


class ObjectManagerCachingRequest(Manager):

    def __init__(self):
        super().__init__()
        if not settings.UNIT_TEST:
            self._queryset_class = QuerySetRequestCache
