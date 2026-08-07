import json
import logging
import os
import re
import time
import traceback
from collections import Counter
from collections.abc import Mapping
from contextlib import ExitStack

from django.db import connection
from django.test import Client, TestCase, override_settings
from django.test.utils import CaptureQueriesContext
from django.urls import reverse
from django.utils.http import urlencode

QUERY_PROFILE_FILE = os.environ.get("VG_QUERY_PROFILE")
QUERY_TRACE_PATTERN = os.environ.get("VG_QUERY_TRACE")  # regex: log stack traces for matching SQL

# Models whose managers (ObjectManagerCachingImmutable/Request) cache lookups in production
# but disable caching under settings.UNIT_TEST - so repeated gets on these tables in a test
# are not real production queries. Used by production_query_count().
PRODUCTION_CACHED_TABLES = (
    "snpdb_genomebuild",
    "genes_genesymbol",
    "flags_flagtype",
    "classification_resolvedvariantinfo",
    "snpdb_allele",
    "snpdb_organization",
    "snpdb_lab",
)


def production_query_count(captured_queries) -> int:
    """ Number of queries that would hit the database in production: excludes test
        savepoints and lookups on tables whose object managers cache in production """
    count = 0
    for query in captured_queries:
        sql = query["sql"]
        if sql.startswith(("SAVEPOINT", "RELEASE SAVEPOINT", "ROLLBACK TO SAVEPOINT")):
            continue
        if sql.startswith("SELECT") and any(f'FROM "{table}"' in sql for table in PRODUCTION_CACHED_TABLES):
            continue
        count += 1
    return count


def _normalize_sql(sql: str) -> str:
    """ Strip literals so queries differing only by parameters group together (N+1 detection) """
    sql = re.sub(r"'[^']*'", "?", sql)
    sql = re.sub(r"\b\d+(\.\d+)?\b", "?", sql)
    sql = re.sub(r"IN \([^)]*\)", "IN (?)", sql)
    return sql


class QueryProfilingClient(Client):
    """ Captures SQL for each GET and appends a JSON line to VG_QUERY_PROFILE """

    def _trace_wrapper(self, path):
        def wrapper(execute, sql, params, many, context):
            if re.search(QUERY_TRACE_PATTERN, sql):
                stack = [line for line in traceback.format_stack()
                         if "/site-packages/" not in line and "unittest_utils" not in line]
                with open(QUERY_PROFILE_FILE + ".trace", "a") as f:
                    f.write(json.dumps({"path": path, "sql": sql[:300], "stack": stack[-12:]}) + "\n")
            return execute(sql, params, many, context)
        return wrapper

    def get(self, path, *args, **kwargs):
        start = time.monotonic()
        with ExitStack() as stack:
            ctx = stack.enter_context(CaptureQueriesContext(connection))
            if QUERY_TRACE_PATTERN:
                stack.enter_context(connection.execute_wrapper(self._trace_wrapper(path)))
            response = super().get(path, *args, **kwargs)
        request_ms = (time.monotonic() - start) * 1000
        real_queries = [q for q in ctx.captured_queries
                        if not q["sql"].startswith(("SAVEPOINT", "RELEASE SAVEPOINT", "ROLLBACK TO SAVEPOINT"))]
        queries = [q["sql"] for q in real_queries]
        duplicates = [{"count": count, "sql": sql}
                      for sql, count in Counter(_normalize_sql(sql) for sql in queries).most_common()
                      if count > 1]
        record = {
            "path": path,
            "status": response.status_code,
            "num_queries": len(queries),
            "request_ms": round(request_ms, 1),
            "sql_ms": round(sum(float(q["time"]) for q in real_queries) * 1000, 1),
            "duplicates": duplicates,
        }
        with open(QUERY_PROFILE_FILE, "a") as f:
            f.write(json.dumps(record) + "\n")
        return response


def _make_test_client() -> Client:
    if QUERY_PROFILE_FILE:
        return QueryProfilingClient()
    return Client()


def prevent_request_warnings(original_function):
    """
    If we need to test for 404s or 405s this decorator can prevent the
    request class from throwing warnings.

    @see https://stackoverflow.com/a/46079090/295724
    """

    def new_function(*args, **kwargs):
        # raise logging level to ERROR
        logger = logging.getLogger('django.request')
        previous_logging_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        # trigger original function that would throw warning
        original_function(*args, **kwargs)

        # lower logging level back to previous
        logger.setLevel(previous_logging_level)

    return new_function


@override_settings(STATICFILES_STORAGE='django.contrib.staticfiles.storage.StaticFilesStorage',
                   VARIANT_ZYGOSITY_GLOBAL_COLLECTION="global",
                   LOG_PARTITION_WARNINGS=False,
                   GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID=None,
                   LIFTOVER_CLASSIFICATIONS=False,
                   ANNOTATION_CACHED_WEB_RESOURCES=[],  # So we don't auto load resources in test
                   CELERY_ALWAYS_EAGER=True)  # Don't launch async tasks
class URLTestCase(TestCase):
    """ Need to override settings as ManifestStaticFilesStorage expects staticfiles.json to exist
        and contain the file asked. @see https://stackoverflow.com/a/51580328/295724 """

    def _test_urls(self, names_and_kwargs, user=None, expected_code_override=None):
        client = _make_test_client()
        if user:
            client.force_login(user)

        for name, kwargs, expected_code in names_and_kwargs:
            with self.subTest(url_name=name):
                expected_code = expected_code_override or expected_code
                kwargs = kwargs.copy()  # As we'll pop next
                get_params = kwargs.pop("GET_PARAMS", {})
                url = reverse(name, kwargs=kwargs)
                if get_params:
                    url += "?" + "&".join([f"{k}={v}" for k, v in get_params.items()])
                # print(f"url: {url}")
                response = client.get(url)
                self.assertEqual(response.status_code, expected_code, msg=f"Url '{url}'")

    def _test_autocomplete_urls(self, names_obj_kwargs, user, in_results):
        client = _make_test_client()
        client.force_login(user)

        for name, obj, get_kwargs in names_obj_kwargs:
            url = reverse(name)
            if get_kwargs:
                url += "?" + urlencode(get_kwargs)
            response = client.get(url)
            self.assertEqual(response.status_code, 200, msg=f"{url} OK")
            if response.status_code == 200:
                data = json.loads(response.content)
                found = False
                for result in data["results"]:
                    if result["id"] == str(obj.pk):
                        found = True
                        break

                #print(f"Url '{url} obj pk={obj.pk} in results={in_results}'")
                self.assertEqual(in_results, found, msg=f"Url '{url} obj pk={obj.pk} in results={in_results}'")

    def _test_datatable_urls(self, names_and_kwargs, user=None, expected_code_override=None):
        """ This calls both the URL and "?dataTableDefinition=1" """

        self._test_urls(names_and_kwargs, user=user, expected_code_override=expected_code_override)

        datatable_definition_names_and_kwargs = []
        for name, kwargs, expected_code in names_and_kwargs:
            new_kwargs = kwargs.copy()
            get_params = new_kwargs.get("GET_PARAMS", {})
            get_params["dataTableDefinition"] = 1
            new_kwargs["GET_PARAMS"] = get_params
            datatable_definition_names_and_kwargs.append((name, new_kwargs, expected_code))

        self._test_urls(datatable_definition_names_and_kwargs, user=user, expected_code_override=expected_code_override)

    def _test_datatables_grid_urls_contains_objs(self, names_obj, user, in_results):
        client = _make_test_client()
        client.force_login(user)

        for name, kwargs, obj in names_obj:
            kwargs = kwargs.copy()  # In case client shared them
            url = reverse(name, kwargs=kwargs)
            definition_url = url + "?dataTableDefinition=1"
            response = client.get(definition_url)
            if in_results:
                response.json()  # Just need to verify that
                self.assertEqual(200, response.status_code)

            response = client.get(url)
            if in_results:
                self.assertEqual(200, response.status_code)
            elif response.status_code == 403:
                continue  # No need to check

            data = json.loads(response.content)
            if obj is not None:
                found = False
                key = "id"
                if isinstance(obj, tuple):
                    key, obj = obj

                for row in data["data"]:
                    if value := row.get("id"):
                        # If what we retrieved was JSON, pull out the key - otherwise just use the value
                        if isinstance(value, Mapping):
                            value = value.get(key)
                    else:
                        value = row.get(key)

                    if value == obj.pk:
                        found = True
                        break

                # print(f"{data=}")
                # print(f"Url '{url} obj pk={obj.pk} in results={in_results}'")
                self.assertEqual(in_results, found, msg=f"Url '{url} obj {key}={obj.pk} in results={in_results}'")

    def _test_jqgrid_urls_contains_objs(self, names_obj, user, in_results):
        """ Also allow 403 if not expecting results
            TODO: Load grid properly call URL with sidx params etc, currently get UnorderedObjectListWarning """
        client = _make_test_client()
        client.force_login(user)

        for name, kwargs, obj in names_obj:
            kwargs = kwargs.copy()  # In case client shared them
            kwargs["op"] = "config"
            config_url = reverse(name, kwargs=kwargs)
            response = client.get(config_url)
            if in_results:
                self.assertEqual(200, response.status_code)
            try:
                config_data = response.json()
            except ValueError:  # Not JSON
                config_data = {}
            kwargs["op"] = "handler"  # To show grid
            url = reverse(name, kwargs=kwargs)
            # Add sidx so we don't get pager order warning
            if sortname := config_data.get("sortname"):
                url += f"?sidx={sortname}"
            response = client.get(url)
            if in_results:
                self.assertEqual(200, response.status_code)
            elif response.status_code == 403:
                continue  # No need to check

            data = json.loads(response.content)
            if obj is not None:
                found = False
                key = "id"
                if isinstance(obj, tuple):
                    key, obj = obj

                for row in data["rows"]:
                    if row.get(key) == obj.pk:
                        found = True
                        break

                # print(f"Url '{url} obj pk={obj.pk} in results={in_results}'")
                self.assertEqual(in_results, found, msg=f"Url '{url} obj {key}={obj.pk} in results={in_results}'")
