from datetime import timedelta

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from snpdb.models import UserSettings
from snpdb.models.models_user_settings import UserSettingsOverride


class VariantTagStaleDateTest(TestCase):
    """ variant_tag_stale_date derives the shared "stale before this" cutoff from
        variant_tag_stale_days (None = staleness off) - #1433 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create(username="tag_stale_user")
        cls.user_override = UserSettingsOverride.objects.get_or_create(user=cls.user)[0]

    def test_none_when_unset(self):
        user_settings = UserSettings.get_for_user(self.user)
        self.assertIsNone(user_settings.variant_tag_stale_days)
        self.assertIsNone(user_settings.variant_tag_stale_date)

    def test_date_when_set(self):
        self.user_override.variant_tag_stale_days = 730
        self.user_override.save()

        user_settings = UserSettings.get_for_user(self.user)
        self.assertEqual(user_settings.variant_tag_stale_days, 730)
        expected = timezone.now() - timedelta(days=730)
        self.assertAlmostEqual(user_settings.variant_tag_stale_date.timestamp(), expected.timestamp(), delta=60)
