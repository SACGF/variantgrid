import json
import logging

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
                get_params = kwargs.pop("GET_PARAMS", {})
                url = reverse(name, kwargs=kwargs)
                if get_params:
                    url += "?" + "&".join([f"{k}={v}" for k, v in get_params.items()])
                    print(f"url: {url}")
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
            data = json.loads(response.content)
            found = False
            for result in data["results"]:
                if result["id"] == str(obj.pk):
                    found = True
                    break

            #print(f"Url '{url} obj pk={obj.pk} in results={in_results}'")
            self.assertEqual(in_results, found, msg=f"Url '{url} obj pk={obj.pk} in results={in_results}'")

    def _test_grid_list_urls(self, names_obj, user, in_results):
        """ Also allow 403 if not expecting results """
        client = Client()
        client.force_login(user)

        for name, kwargs, obj in names_obj:
            kwargs = kwargs.copy()  # In case client shared them
            kwargs["op"] = "handler"  # To show grid
            url = reverse(name, kwargs=kwargs)
            response = client.get(url)
            if response.status_code == 403 and in_results is False:
                return
            self.assertEqual(200, response.status_code)

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

                #print(f"Url '{url} obj pk={obj.pk} in results={in_results}'")
                self.assertEqual(in_results, found, msg=f"Url '{url} obj {key}={obj.pk} in results={in_results}'")
