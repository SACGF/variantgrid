"""
Redis cache backend that treats an entry it can no longer unpickle as a miss: the key is deleted and the
default returned. An entry pickled before a deploy that renamed, moved or removed its class then costs one
recompute instead of raising on every read until `CACHE_VERSION` is bumped.

Entry point: `UnpickleTolerantRedisCache`, the `BACKEND` of `CACHES` in
`variantgrid/settings/components/default_settings.py`.
"""
import logging
import pickle

from django.core.cache.backends.redis import RedisCache, RedisCacheClient

# What pickle.loads raises when the class an entry names has gone: renamed or removed (AttributeError),
# moved to another module (ModuleNotFoundError, an ImportError), or bytes it can't parse
UNPICKLE_ERRORS = (pickle.UnpicklingError, AttributeError, ImportError, EOFError)

_MISSING = object()


class UnpickleTolerantRedisCacheClient(RedisCacheClient):
    def get(self, key, default):
        try:
            return super().get(key, default)
        except UNPICKLE_ERRORS as e:
            logging.warning("Deleting cache entry '%s' that can't be unpickled: %s", key, e)
            self.delete(key)
            return default

    def get_many(self, keys):
        try:
            return super().get_many(keys)
        except UNPICKLE_ERRORS:
            # One stale entry fails the whole batch - read each so only the bad ones are dropped
            values = {k: self.get(k, _MISSING) for k in keys}
            return {k: v for k, v in values.items() if v is not _MISSING}


class UnpickleTolerantRedisCache(RedisCache):
    def __init__(self, server, params):
        super().__init__(server, params)
        self._class = UnpickleTolerantRedisCacheClient
