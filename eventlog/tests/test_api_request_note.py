from django.test import SimpleTestCase

from eventlog.models import api_request_note
from library.django_utils.change_tracking import collect_changes


class ApiRequestNoteTest(SimpleTestCase):

    def test_note_goes_to_the_active_request_only(self):
        api_request_note("outside a tracked request")  # only logged
        with collect_changes() as collector:
            api_request_note("inside")
        self.assertEqual(["inside"], collector.notes)
