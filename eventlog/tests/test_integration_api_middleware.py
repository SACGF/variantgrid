from django.http import HttpResponse
from django.test import RequestFactory, TestCase, override_settings
from django.urls import reverse

from eventlog.middleware import IntegrationApiMiddleware
from eventlog.models import IntegrationActivity

TRACKED_PATH = reverse('eventlog_detail', kwargs={'pk': 1})
TRACKED_PREFIX = TRACKED_PATH.rsplit("/", 1)[0]
TRACKING = {TRACKED_PREFIX: "Lab pipeline API"}


class IntegrationApiMiddlewareTest(TestCase):

    def setUp(self):
        self.factory = RequestFactory()

    def _call(self, path: str, status_code: int = 200):
        middleware = IntegrationApiMiddleware(lambda request: HttpResponse(status=status_code))
        middleware(self.factory.post(path))

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
