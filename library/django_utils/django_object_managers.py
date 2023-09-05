from dataclasses import dataclass
from typing import Dict, Any, Optional

from django.db.models import Manager, QuerySet, Q
from frozendict import frozendict


@dataclass(frozen=True)
class CachedQuery:
    model: Any
    args: Any
    kwargs: Any


_CACHED_QUERIES: Dict[CachedQuery, Any] = {}
# TODO, could make this a LRU cache so it doesn't grow infinitely for larger models


class QuerySetCaching(QuerySet):

    def __init__(self, model=None, query=None, using=None, hints=None):
        super().__init__(model=model, query=query, using=using, hints=hints)

    def get(self, *args, **kwargs):
        cq: Optional[CachedQuery] = None
        kwargs_fd = None
        if kwargs:
            # turn into a frozen dict so we can hash it
            kwargs_fd = frozendict(kwargs)
        try:
            cq = CachedQuery(model=self.model, args=args, kwargs=kwargs_fd)
        except TypeError:
            # unhashable argument
            pass

        if cq:
            if existing := _CACHED_QUERIES.get(cq):
                return existing
            else:
                result = super().get(*args, **kwargs)
                _CACHED_QUERIES[cq] = result
                return result
        else:
            return super().get(*args, **kwargs)


class CachingObjectManager(Manager):
    """
    Best used as both the default objects and as Meta._base_manager_name = 'objects'
    Will cache results for queries - so only do this on models that are read from but never written to.

    We could invalidate the cache when save() is called, but that only works if the save is done in the only instance
    that has already cached values.

    Note that currently it only caches the results of get(...)
    """

    def __init__(self):
        super().__init__()
        self._queryset_class = QuerySetCaching
