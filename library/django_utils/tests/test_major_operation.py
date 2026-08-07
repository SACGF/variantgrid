"""
Tests for the major_operation per-user concurrency framework.

@see variantgrid_private #1502
"""
from contextlib import ExitStack
from unittest import mock

from django.core.cache import cache
from django.test import TestCase, override_settings

from library.django_utils.major_operation import major_operation, TooManyMajorOperationsError

LOCMEM_CACHE = {
    "default": {"BACKEND": "django.core.cache.backends.locmem.LocMemCache"},
}


# Patch out the Postgres-only statement_timeout so these run against any/no DB connection
@mock.patch("library.django_utils.major_operation._statement_timeout")
@override_settings(CACHES=LOCMEM_CACHE, MAJOR_OPERATION_LIMITS_ENABLED=True,
                   MAJOR_OPERATION_MAX_CONCURRENT_PER_USER=3, MAJOR_OPERATION_SLOT_EXPIRE_SECONDS=600)
class MajorOperationTests(TestCase):
    USER = "test_user"

    def setUp(self):
        cache.clear()

    def test_under_limit_runs_and_releases(self, _mock_timeout):
        with major_operation(self.USER, "grid"):
            pass
        # Slot released - we can run the limit again afterwards
        for _ in range(3):
            with major_operation(self.USER, "grid"):
                pass

    def test_exceeding_limit_raises(self, _mock_timeout):
        with ExitStack() as stack:
            for _ in range(3):  # fill all 3 slots
                stack.enter_context(major_operation(self.USER, "grid"))
            with self.assertRaises(TooManyMajorOperationsError):
                stack.enter_context(major_operation(self.USER, "grid"))

    def test_slot_freed_lets_new_operation_run(self, _mock_timeout):
        with ExitStack() as stack:
            for _ in range(3):
                stack.enter_context(major_operation(self.USER, "grid"))
            # Free one slot
            stack.close()
        # Now a fresh operation succeeds
        with major_operation(self.USER, "grid"):
            pass

    def test_failed_operation_releases_slot(self, _mock_timeout):
        with self.assertRaises(ValueError):
            with major_operation(self.USER, "grid"):
                raise ValueError("boom")
        # The slot from the failed op must have been released - 3 more should fit
        with ExitStack() as stack:
            for _ in range(3):
                stack.enter_context(major_operation(self.USER, "grid"))

    def test_per_user_isolation(self, _mock_timeout):
        with ExitStack() as stack:
            for _ in range(3):  # fill user A's slots
                stack.enter_context(major_operation("user_a", "grid"))
            # user B is unaffected
            with major_operation("user_b", "grid"):
                pass

    @override_settings(MAJOR_OPERATION_LIMITS_ENABLED=False)
    def test_disabled_does_not_limit(self, _mock_timeout):
        with ExitStack() as stack:
            for _ in range(10):
                stack.enter_context(major_operation(self.USER, "grid"))
