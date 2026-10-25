import json

from django.contrib.auth.models import User
from django.http import HttpResponse
from django.test import RequestFactory, TestCase, override_settings
from django.urls import reverse

from eventlog.middleware import IntegrationApiMiddleware
from eventlog.models import ApiRequestLog, IntegrationActivity

TRACKED_PATH = reverse('eventlog_detail', kwargs={'pk': 1})
TRACKED_PREFIX = TRACKED_PATH.rsplit("/", 1)[0]
TRACKING = {TRACKED_PREFIX: "Lab pipeline API"}


class IntegrationApiMiddlewareTest(TestCase):

    def setUp(self):
        self.factory = RequestFactory()

    def _call(self, path: str, status_code: int = 200, request=None, content=b""):
        middleware = IntegrationApiMiddleware(lambda _request: HttpResponse(content, status=status_code))
        middleware(request or self.factory.post(path))

    @override_settings(INTEGRATION_API_TRACKING=TRACKING)
    def test_records_against_a_configured_prefix(self):
        self._call(TRACKED_PATH)
        activity = IntegrationActivity.objects.get(key="api:eventlog_detail")
        self.assertEqual("Lab pipeline API (eventlog_detail)", activity.name)
        self.assertIsNotNone(activity.last_success)

    @override_settings(INTEGRATION_API_TRACKING=TRACKING)
    def test_ignores_an_unconfigured_prefix(self):
        self._call(reverse('eventlog'))
        self.assertFalse(IntegrationActivity.objects.exists())
        self.assertFalse(ApiRequestLog.objects.exists())

    @override_settings(INTEGRATION_API_TRACKING={})
    def test_records_nothing_when_unconfigured(self):
        self._call(TRACKED_PATH)
        self.assertFalse(IntegrationActivity.objects.exists())

    @override_settings(INTEGRATION_API_TRACKING=TRACKING)
    def test_records_an_error_for_a_4xx(self):
        self._call(TRACKED_PATH, status_code=403)
        activity = IntegrationActivity.objects.get(key="api:eventlog_detail")
        self.assertIsNone(activity.last_success)
        self.assertIsNotNone(activity.last_error)
        self.assertIn("403", activity.last_error_message)

    @override_settings(INTEGRATION_API_TRACKING=TRACKING)
    def test_logs_the_call_with_secrets_redacted(self):
        user = User.objects.create_user("api_log_pipeline")
        data = {"sequencing_run": "run_1", "auth": {"api_token": "abc"}, "samples": [{"password": "x", "name": "s1"}]}
        request = self.factory.post(TRACKED_PATH, data=json.dumps(data), content_type="application/json")
        request.user = user
        self._call(TRACKED_PATH, request=request)

        api_request = ApiRequestLog.objects.get()
        self.assertEqual(user, api_request.user)
        self.assertEqual(200, api_request.status_code)
        self.assertEqual({"pk": 1}, api_request.view_kwargs)
        logged = json.loads(api_request.request_body)
        self.assertEqual("run_1", logged["sequencing_run"])
        self.assertEqual("[redacted]", logged["auth"]["api_token"])
        self.assertEqual({"password": "[redacted]", "name": "s1"}, logged["samples"][0])
        self.assertIsNone(api_request.response_body)

    @override_settings(INTEGRATION_API_TRACKING=TRACKING)
    def test_form_body_and_query_string_redacted(self):
        request = self.factory.post(f"{TRACKED_PATH}?api_key=abc&run=1", data="password=x&name=s1",
                                    content_type="application/x-www-form-urlencoded")
        self._call(TRACKED_PATH, request=request)
        api_request = ApiRequestLog.objects.get()
        self.assertEqual("api_key=%5Bredacted%5D&run=1", api_request.query_string)
        self.assertEqual("password=%5Bredacted%5D&name=s1", api_request.request_body)

    @override_settings(INTEGRATION_API_TRACKING=TRACKING)
    def test_multipart_body_is_noted_not_read(self):
        request = self.factory.post(TRACKED_PATH, data={"name": "value"})
        self._call(TRACKED_PATH, request=request)
        self.assertFalse(hasattr(request, "_body"))
        self.assertTrue(ApiRequestLog.objects.get().request_body.startswith("[multipart/form-data, "))

    @override_settings(INTEGRATION_API_TRACKING=TRACKING, API_REQUEST_LOG_MAX_BODY=10)
    def test_keeps_a_trimmed_response_body_for_an_error(self):
        self._call(TRACKED_PATH, status_code=400, content=b'{"sequencing_run": ["Not found"]}')
        api_request = ApiRequestLog.objects.get()
        self.assertTrue(api_request.response_body.startswith('{"sequenci\n[trimmed from 33 characters]'))
