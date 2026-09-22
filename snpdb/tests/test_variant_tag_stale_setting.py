from datetime import timedelta
from importlib import import_module

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from snpdb.models import TagConfigCollection, UserSettings
from snpdb.models.models_user_settings import UserSettingsOverride

# The migration module name starts with a digit, so it can't be imported with a plain import statement
copy_stale_days_to_collections = import_module(
    "snpdb.migrations.0273_rename_tagcolor_tagconfig").copy_stale_days_to_collections


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


class FakeCollection:
    def __init__(self, pk):
        self.pk = pk
        self.variant_tag_stale_days = None

    def save(self):
        pass


class FakeOverride:
    def __init__(self, collection, variant_tag_stale_days):
        self.tag_config = collection
        self.variant_tag_stale_days = variant_tag_stale_days


class FakeManager:
    def __init__(self, overrides):
        self.overrides = overrides

    def filter(self, **_kwargs):
        return [o for o in self.overrides if o.tag_config and o.variant_tag_stale_days is not None]


class FakeModel:
    def __init__(self, overrides):
        self.objects = FakeManager(overrides)


class FakeApps:
    """ Stands in for the migration's historical models - only objects.filter() is used """

    def __init__(self, overrides_by_model):
        self.overrides_by_model = overrides_by_model

    def get_model(self, _app_label, model_name):
        return FakeModel(self.overrides_by_model.get(model_name, []))


class CopyStaleDaysMigrationTest(TestCase):
    """ 0273's data copy: the most specific override's window is what lands on the collection,
        and one user-level override per collection wins """

    def test_user_override_beats_lab(self):
        collection = FakeCollection(1)
        apps = FakeApps({
            "LabUserSettingsOverride": [FakeOverride(collection, 365)],
            "UserSettingsOverride": [FakeOverride(collection, 730)],
        })
        copy_stale_days_to_collections(apps)
        self.assertEqual(collection.variant_tag_stale_days, 730)

    def test_first_user_override_keeps_shared_collection(self):
        collection = FakeCollection(1)
        apps = FakeApps({
            "UserSettingsOverride": [FakeOverride(collection, 730), FakeOverride(collection, 180)],
        })
        copy_stale_days_to_collections(apps)
        self.assertEqual(collection.variant_tag_stale_days, 730)
