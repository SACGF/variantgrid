from django.contrib.auth.models import User
from django.test import TestCase

from flags.models.enums import FlagStatus
from flags.models.models import (
    FlagCollection,
    FlagResolution,
    FlagType,
    FlagTypeContext,
    FlagTypeResolution,
)


class FlagPermissionTest(TestCase):
    """ A FlagCollection with no source object gives NO_PERM to everyone but superusers (ADMIN) and admin_bot """

    @classmethod
    def setUpTestData(cls):
        context = FlagTypeContext.objects.create(id='test_flag_permissions', label='Test Flag Permissions')
        cls.flag_type = FlagType.objects.create(id='test_flag_permissions_users', context=context, label='Test',
                                                description='Test', permission='O', raise_permission='U')
        open_resolution, _ = FlagResolution.objects.get_or_create(id='open', defaults={
            'label': 'In Progress', 'description': 'This is ongoing', 'status': FlagStatus.OPEN})
        FlagTypeResolution.objects.create(flag_type=cls.flag_type, resolution=open_resolution)
        cls.collection = FlagCollection.objects.create(context=context)
        cls.user = User.objects.create_user(username='test_flag_permissions_user')
        cls.superuser = User.objects.create_superuser(username='test_flag_permissions_superuser')

    def test_add_flag_permission_check(self):
        with self.assertRaises(PermissionError):
            self.collection.add_flag(self.flag_type, user=self.user, permission_check=True)
        flag = self.collection.add_flag(self.flag_type, user=self.superuser, permission_check=True)
        self.assertEqual(flag.user, self.superuser)

    def test_flags_hidden_from_no_perm_user(self):
        flag = self.collection.add_flag(self.flag_type)
        self.assertFalse(self.collection.flags(user=self.user).exists())
        self.assertEqual(list(self.collection.flags(user=self.superuser)), [flag])
