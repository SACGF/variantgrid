from unittest.mock import patch

import redis
from django.core.cache import cache
from django.test import SimpleTestCase

from upload.tasks.vcf.unknown_variants_task import InsertUnknownVariantsTask


@patch("upload.tasks.vcf.unknown_variants_task.time.sleep")
class TestInsertUnknownVariantsLock(SimpleTestCase):
    """ #2189 - a held lock or a Redis outage must wait, never skip the insert """

    def tearDown(self):
        cache.delete(InsertUnknownVariantsTask.LOCK_ID)

    def test_held_lock_raises_rather_than_skipping(self, _sleep):
        cache.add(InsertUnknownVariantsTask.LOCK_ID, "true", 60)
        with patch.object(InsertUnknownVariantsTask, "LOCK_WAIT", -1):
            with self.assertRaises(TimeoutError):
                InsertUnknownVariantsTask._acquire_lock()

    def test_redis_loading_retries_then_acquires(self, sleep):
        real_add = cache.add
        responses = [redis.exceptions.BusyLoadingError("loading"), None]

        def flaky_add(*args, **kwargs):
            if responses and (error := responses.pop(0)):
                raise error
            return real_add(*args, **kwargs)

        with patch.object(cache, "add", side_effect=flaky_add):
            InsertUnknownVariantsTask._acquire_lock()
        self.assertEqual(sleep.call_count, 1)
        self.assertIsNotNone(cache.get(InsertUnknownVariantsTask.LOCK_ID))
