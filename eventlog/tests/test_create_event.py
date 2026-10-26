""" create_event view: only records the events the browser JS sends """
from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from eventlog.models import Event


class CreateEventTest(TestCase):

    def setUp(self):
        self.client.force_login(User.objects.create_user("create_event_user"))

    def test_only_client_events_recorded(self):
        js_error = {"app_name": "analysis", "event_name": "javascript_error", "severity": "E", "details": "x"}
        for data, status_code in [
            (js_error, 200),
            ({**js_error, "app_name": "upload", "event_name": "import_failed"}, 400),
            ({**js_error, "severity": "I"}, 400),
        ]:
            with self.subTest(data=data):
                self.assertEqual(status_code, self.client.post(reverse("create_event"), data).status_code)
        self.assertEqual(1, Event.objects.exclude(name="login").count())
