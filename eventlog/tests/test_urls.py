import unittest

from django.contrib.auth.models import User

from eventlog.models import create_event
from library.django_utils.unittest_utils import URLTestCase


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        owner_username = f"test_user_{__file__}_owner"
        admin_username = f"test_user_{__file__}_admin"
        cls.user = User.objects.get_or_create(username=owner_username)[0]
        cls.admin_user = User.objects.create_superuser(admin_username)
        cls.event = create_event(user=cls.user, name="Test Event")

    def testDataGridUrls(self):
        DATATABLE_GRID_LIST_URLS = [
            ("event_log_datatable", {}, 200),
        ]
        self._test_datatable_urls(DATATABLE_GRID_LIST_URLS, self.user)

    def testApiRequestLogUrls(self):
        self._test_datatable_urls([("api_request_log_datatable", {}, 200)], self.admin_user)
        self._test_datatable_urls([("api_request_log_datatable", {}, 403)], self.user)
        object_kwargs = {"app_label": "seqauto", "model_name": "sequencingrun", "object_pk": "RUN_1"}
        self._test_datatable_urls([("api_request_log_object_datatable", object_kwargs, 200)], self.admin_user)
        self._test_datatable_urls([("api_request_log_object_datatable", object_kwargs, 403)], self.user)

if __name__ == "__main__":
    unittest.main()
