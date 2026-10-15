"""
The predicate the Redis OOM manual migration and deployment_check share: when a box needs the
OOMScoreAdjust drop-in for redis (#2187).
"""
from unittest.mock import patch

import redis
from django.test import SimpleTestCase

from variantgrid.deployment_validation.memory_check import redis_unprotected_from_oom_killer

OOM_SCORE_ADJ = "variantgrid.deployment_validation.memory_check._redis_oom_score_adj"


class RedisUnprotectedFromOomKillerTest(SimpleTestCase):

    def test_default_score_needs_the_drop_in(self):
        with patch(OOM_SCORE_ADJ, return_value=0):
            self.assertTrue(redis_unprotected_from_oom_killer())

    def test_lowered_score_is_protected(self):
        with patch(OOM_SCORE_ADJ, return_value=-900):
            self.assertFalse(redis_unprotected_from_oom_killer())

    def test_remote_redis_has_nothing_to_install(self):
        with patch(OOM_SCORE_ADJ, return_value=None):
            self.assertFalse(redis_unprotected_from_oom_killer())

    def test_unreachable_redis_does_not_fail_migrate(self):
        with patch(OOM_SCORE_ADJ, side_effect=redis.ConnectionError("refused")):
            self.assertFalse(redis_unprotected_from_oom_killer())
