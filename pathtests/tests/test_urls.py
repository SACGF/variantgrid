from django.contrib.auth.models import User

from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from library.guardian_utils import assign_permission_to_user_and_groups
from pathtests.models import Case, PathologyTestOrder
from patients.models import Patient


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user_owner = User.objects.get_or_create(username='pathtests_user')[0]
        cls.user_non_owner = User.objects.get_or_create(username='pathtests_other_user')[0]
        patient = Patient.objects.create(first_name="Path", last_name="Tests")
        assign_permission_to_user_and_groups(cls.user_owner, patient)
        cls.case = Case.objects.create(name="pathtests_case", patient=patient)
        cls.pathology_test_order = PathologyTestOrder.objects.create(case=cls.case)

        cls.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS = [
            ("view_case", {"pk": cls.case.pk}, 200),
            ("view_pathology_test_order", {"pk": cls.pathology_test_order.pk}, 200),
        ]
        cls.PRIVATE_DATATABLES_GRID_LIST_URLS = [
            ("cases_datatable", {}, ("text", cls.case)),
            ("pathology_test_orders_datatable", {}, ("text", cls.pathology_test_order)),
        ]
        cls.PRIVATE_AUTOCOMPLETE_URLS = [
            ('case_autocomplete', cls.case, {"q": cls.case.name}),
        ]

    def testPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_owner)

    @prevent_request_warnings
    def testNoPermission(self):
        self._test_urls(self.PRIVATE_OBJECT_URL_NAMES_AND_KWARGS, self.user_non_owner, expected_code_override=403)

    def testDatatableUrls(self):
        DATATABLE_URLS = [
            ("pathology_test_orders_datatable", {}, 200),
            ("cases_datatable", {}, 200),
            ("pathology_tests_datatable", {}, 200),
        ]
        self._test_datatable_urls(DATATABLE_URLS, self.user_non_owner)

    def testGridListPermission(self):
        self._test_datatables_grid_urls_contains_objs(self.PRIVATE_DATATABLES_GRID_LIST_URLS, self.user_owner, True)

    def testGridListNoPermission(self):
        self._test_datatables_grid_urls_contains_objs(self.PRIVATE_DATATABLES_GRID_LIST_URLS, self.user_non_owner, False)

    def testAutocompletePermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_owner, True)

    def testAutocompleteNoPermission(self):
        self._test_autocomplete_urls(self.PRIVATE_AUTOCOMPLETE_URLS, self.user_non_owner, False)
