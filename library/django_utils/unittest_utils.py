import json
import logging
from collections.abc import Mapping

from django.test import Client, TestCase, override_settings
from django.urls import reverse
from django.utils.http import urlencode


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
        client = Client()
        if user:
            client.force_login(user)

        for name, kwargs, expected_code in names_and_kwargs:
            with self.subTest(url_name=name):
                expected_code = expected_code_override or expected_code
                kwargs = kwargs.copy() # As we'll pop next
                get_params = kwargs.pop("GET_PARAMS", {})
                url = reverse(name, kwargs=kwargs)
                if get_params:
                    url += "?" + "&".join([f"{k}={v}" for k, v in get_params.items()])
                # print(f"url: {url}")
                response = client.get(url)
                self.assertEqual(response.status_code, expected_code, msg=f"Url '{url}'")

    def _test_autocomplete_urls(self, names_obj_kwargs, user, in_results):
        client = Client()
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
        client = Client()
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
                        # If the "id" object is JSON, pull out the key - otherwise just use the value
                        if isinstance(value, Mapping):
                            value = value.get(key)
                        if value == obj.pk:
                            found = True
                            break

                # print(f"{data=}")
                # print(f"Url '{url} obj pk={obj.pk} in results={in_results}'")
                self.assertEqual(in_results, found, msg=f"Url '{url} obj {key}={obj.pk} in results={in_results}'")


    def _test_jqgrid_urls_contains_objs(self, names_obj, user, in_results):
        """ Also allow 403 if not expecting results
            TODO: Load grid properly call URL with sidx params etc, currently get UnorderedObjectListWarning """
        client = Client()
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
