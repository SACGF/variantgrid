import unittest

from django.contrib.auth.models import User
from manual.models.manual_migration_models import ManualMigrationTask, ManualMigrationAttempt
from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings


class Test(URLTestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        owner_username = f"test_user_{__file__}_owner"
        admin_username = f"test_user_{__file__}_admin"
        cls.user = User.objects.get_or_create(username=owner_username)[0]
        cls.admin_user = User.objects.create_superuser(admin_username)

        # need to create a manual migration task (and describe it) and manual migration event
        cls.manual_migration_task = ManualMigrationTask.objects.create(id="a*test*migartion")
        cls.manual_migration_event = ManualMigrationAttempt.objects.create(task=cls.manual_migration_task)

        cls.PRIVATE_DATATABLES_GRID_LIST_URLS = [
            ("manual_migrations_datatable", {}, cls.manual_migration_event)
        ]

    @prevent_request_warnings
    def testDataTablesGridListNoPermission(self):
        self._test_datatables_grid_urls_contains_objs(self.PRIVATE_DATATABLES_GRID_LIST_URLS, self.user, True)

if __name__ == "__main__":
    unittest.main()
