"""
"Major operations" are expensive requests (eg loading an analysis node grid over millions of
variants) that can hammer the database. To stop a single user - whether impatient or malicious -
from launching a flood of these concurrently and thrashing Postgres, we:

  * cap how many major operations a user can have running at once (per-user concurrency limit)
  * lower the DB ``statement_timeout`` for the duration of the operation so runaway queries die sooner

Wrap expensive work in the ``major_operation`` context manager. It raises
``TooManyMajorOperationsError`` when the per-user limit is exceeded - callers decide how to respond
(eg redirect-and-retry for grids, or HTTP 503 for APIs).

See variantgrid_private #1502.
"""
import logging
from contextlib import contextmanager

from django.conf import settings
from django.core.cache import cache
from django.db import connection

from library.utils.hash_utils import sha256sum_str


class TooManyMajorOperationsError(Exception):
    """ Raised when a user already has the maximum number of concurrent major operations running """

    def __init__(self, user, operation_name: str, limit: int):
        self.user = user
        self.operation_name = operation_name
        self.limit = limit
        super().__init__(f"User '{user}' exceeded the limit of {limit} concurrent major operations "
                         f"(attempted '{operation_name}')")


def _major_operation_count_key(user) -> str:
    # Single global counter per-user, shared across all major operation types
    return sha256sum_str(f"major_operation_count_{user}")


@contextmanager
def _statement_timeout(seconds: int):
    """ Temporarily lower the Postgres statement_timeout on the current connection.
        Reset on exit as connections are reused (CONN_MAX_AGE). """
    if connection.vendor != 'postgresql':
        yield
        return

    with connection.cursor() as cursor:
        cursor.execute("SET statement_timeout TO %s;", [seconds * 1000])
    try:
        yield
    finally:
        # Restore the global backstop value applied to all connections (see variantgrid/wsgi.py)
        default_ms = settings.DATABASE_STATEMENT_TIMEOUT_SECONDS * 1000
        with connection.cursor() as cursor:
            cursor.execute("SET statement_timeout TO %s;", [default_ms])


@contextmanager
def major_operation(user, operation_name: str = "major_operation"):
    """ Limit concurrent expensive operations per-user and bound their DB query time.

        Raises TooManyMajorOperationsError if the per-user concurrency limit is exceeded.
        operation_name is used for logging/metrics only - the limit is shared across all types. """
    if not settings.MAJOR_OPERATION_LIMITS_ENABLED:
        yield
        return

    limit = settings.MAJOR_OPERATION_MAX_CONCURRENT_PER_USER
    expire = settings.MAJOR_OPERATION_SLOT_EXPIRE_SECONDS
    count_key = _major_operation_count_key(user)

    # Ensure the key exists (Redis incr raises on a missing key) then atomically claim a slot
    cache.add(count_key, 0, expire)
    try:
        current = cache.incr(count_key)
    except ValueError:
        # Key expired between add and incr - reset it
        cache.set(count_key, 1, expire)
        current = 1

    if current > limit:
        cache.decr(count_key)
        logging.warning("User '%s' hit major operation limit of %d (attempted '%s')",
                        user, limit, operation_name)
        raise TooManyMajorOperationsError(user, operation_name, limit)

    try:
        with _statement_timeout(settings.MAJOR_OPERATION_STATEMENT_TIMEOUT_SECONDS):
            yield
    finally:
        cache.decr(count_key)
