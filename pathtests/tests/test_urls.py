from django.contrib.auth.models import User

from library.django_utils.unittest_utils import URLTestCase
from pathtests.models import Case
from patients.models import Patient


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='pathtests_user')[0]
        patient = Patient.objects.create(first_name="Path", last_name="Tests")
        cls.case = Case.objects.create(name="pathtests_case", patient=patient)

    def testDatatableUrls(self):
        DATATABLE_URLS = [
            ("pathology_test_orders_datatable", {}, 200),
            ("cases_datatable", {}, 200),
            ("pathology_tests_datatable", {}, 200),
        ]
        self._test_datatable_urls(DATATABLE_URLS, self.user)

    def testAutocompleteUrls(self):
        AUTOCOMPLETE_URLS = [
            ('case_autocomplete', self.case, {"q": self.case.name}),
        ]
        self._test_autocomplete_urls(AUTOCOMPLETE_URLS, self.user, True)
