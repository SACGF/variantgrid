from datetime import timedelta

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from snpdb.models import TagConfigCollection, UserSettings
from snpdb.models.models_user_settings import UserSettingsOverride


class VariantTagStaleDateTest(TestCase):
    """ variant_tag_stale_date derives the shared "stale before this" cutoff from the tag config
        collection's variant_tag_stale_days (no collection / blank = staleness off) - #1433, #1892 """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create(username="tag_stale_user")
        cls.user_override = UserSettingsOverride.objects.get_or_create(user=cls.user)[0]

    def test_none_when_collection_has_no_days(self):
        collection = TagConfigCollection.objects.create(name="no staleness", user=self.user)
        self.user_override.tag_config = collection
        self.user_override.save()

        user_settings = UserSettings.get_for_user(self.user)
        self.assertIsNone(user_settings.variant_tag_stale_date)

    def test_date_from_collection(self):
        collection = TagConfigCollection.objects.create(name="2 years", user=self.user,
                                                        variant_tag_stale_days=730)
        self.user_override.tag_config = collection
        self.user_override.save()

        user_settings = UserSettings.get_for_user(self.user)
        expected = timezone.now() - timedelta(days=730)
        self.assertAlmostEqual(user_settings.variant_tag_stale_date.timestamp(), expected.timestamp(), delta=60)
