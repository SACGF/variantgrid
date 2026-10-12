import pickle
from unittest import mock

from django.test import SimpleTestCase

from library.django_utils.redis_cache import UnpickleTolerantRedisCacheClient


class _Renamed:
    pass


class _FakeRedis:
    def __init__(self, data: dict):
        self.data = data

    def get(self, key):
        return self.data.get(key)

    def mget(self, keys):
        return [self.data.get(k) for k in keys]

    def delete(self, key):
        return self.data.pop(key, None) is not None


class UnpickleTolerantRedisCacheTest(SimpleTestCase):
    def setUp(self):
        # As if _Renamed was renamed to _Missing after the entry was written (same length keeps the pickle valid)
        stale = pickle.dumps(_Renamed()).replace(b"_Renamed", b"_Missing")
        self.redis = _FakeRedis({"good": pickle.dumps("value"), "stale": stale})
        self.client = UnpickleTolerantRedisCacheClient(["redis://unused"])
        mock.patch.object(self.client, "get_client", return_value=self.redis).start()
        self.addCleanup(mock.patch.stopall)

    def test_get_stale_entry_is_a_miss_and_deleted(self):
        self.assertEqual(self.client.get("stale", "default"), "default")
        self.assertNotIn("stale", self.redis.data)
        self.assertEqual(self.client.get("good", None), "value")

    def test_get_many_drops_only_stale_entries(self):
        self.assertEqual(self.client.get_many(["good", "stale", "absent"]), {"good": "value"})
        self.assertNotIn("stale", self.redis.data)
