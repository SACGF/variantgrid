import unittest

from django.contrib.auth.models import User
from django.test.testcases import TestCase

from eventlog.models import create_event


class EventLogModelTestCase(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        owner_username = f"test_user_{__file__}_owner"
        non_owner_username = f"test_user_{__file__}_non_owner"
        admin_username = f"test_user_{__file__}_admin"

        cls.owner_user = User.objects.get_or_create(username=owner_username)[0]
        cls.non_owner_user = User.objects.get_or_create(username=non_owner_username)[0]
        cls.admin_user = User.objects.create_superuser(admin_username)

    def test_permissions(self):
        event = create_event(user=self.owner_user, name="Test Event")
        self.assertTrue(event.can_write(self.owner_user))
        self.assertFalse(event.can_write(self.non_owner_user))

        event = create_event(user=None, name="No user Test Event")
        self.assertTrue(event.can_write(self.admin_user))
        self.assertFalse(event.can_write(self.non_owner_user))



if __name__ == '__main__':
    unittest.main()
