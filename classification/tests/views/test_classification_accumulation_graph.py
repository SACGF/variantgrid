from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse


class AccumulationReportFrequencyTestCase(TestCase):
    """ The report's time step comes from the request, so it is bounded before the report steps through it """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.create_user(username="accumulation_report_user")

    def _get(self, frequency: str):
        self.client.force_login(self.user)
        return self.client.get(reverse("classification_accumulation_data"), {"frequency": frequency})

    def test_frequency_out_of_range_refused(self):
        for frequency in ["-1", "0", "367", "weekly"]:
            with self.subTest(frequency=frequency):
                self.assertEqual(self._get(frequency).status_code, 400)

    def test_frequency_in_range_allowed(self):
        self.assertEqual(self._get("1").status_code, 200)
