from django.contrib.auth.models import User

from library.django_utils.unittest_utils import URLTestCase


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='pathtests_user')[0]

    def testDatatableUrls(self):
        DATATABLE_URLS = [
            ("pathology_test_orders_datatable", {}, 200),
            ("cases_datatable", {}, 200),
            ("pathology_tests_datatable", {}, 200),
        ]
        self._test_datatable_urls(DATATABLE_URLS, self.user)
