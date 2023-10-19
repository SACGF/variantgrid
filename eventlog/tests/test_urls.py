import unittest

from django.contrib.auth.models import User
from eventlog.models import create_event
from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings

class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        owner_username = f"test_user_{__file__}_owner"
        admin_username = f"test_user_{__file__}_admin"
        cls.user = User.objects.get_or_create(username=owner_username)[0]
        cls.admin_user = User.objects.create_superuser(admin_username)
        cls.event = create_event(user=cls.user, name="Test Event")

        cls.PRIVATE_DATATABLES_GRID_LIST_URLS = [
            ("event_log_datatable", {}, cls.event)
        ]

    def testDataGridUrls(self):
        DATATABLE_GRID_LIST_URLS = [
            ("event_log_datatable", {}, 200),
        ]
        self._test_datatable_urls(DATATABLE_GRID_LIST_URLS, self.user)

    @prevent_request_warnings
    def testDataTablesGridListNoPermission(self):
        self._test_datatables_grid_urls_contains_objs(self.PRIVATE_DATATABLES_GRID_LIST_URLS, self.user, True)

if __name__ == "__main__":
    unittest.main()
