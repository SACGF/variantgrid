from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models import UserSettings, Lab, Organization, Country
from snpdb.models.models_user_settings import GlobalSettings, OrganizationUserSettingsOverride, \
    LabUserSettingsOverride, UserSettingsOverride


class NodeGridAutoLoadSettingTest(TestCase):
    """ node_grid_auto_load_max_variants cascades Global -> Org -> Lab -> User,
        where later non-null levels override earlier ones. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.org = Organization.objects.create(name="InstX", group_name="instx")
        country = Country.objects.get_or_create(name="CountryA")[0]
        cls.lab = Lab.objects.create(name="Labby", organization=cls.org, city="CityA",
                                     country=country, group_name="instx/labby")
        cls.user = User.objects.create(username="node_grid_user")

        cls.global_settings = GlobalSettings.objects.get()
        cls.org_override = OrganizationUserSettingsOverride.objects.get_or_create(organization=cls.org)[0]
        cls.lab_override = LabUserSettingsOverride.objects.get_or_create(lab=cls.lab)[0]
        cls.user_override = UserSettingsOverride.objects.get_or_create(user=cls.user)[0]
        # User belongs to the lab and uses it by default, so the full cascade applies
        cls.user.groups.add(cls.lab.group)
        cls.user_override.default_lab = cls.lab
        cls.user_override.save()

    def _resolved(self):
        return UserSettings.get_for_user(self.user).node_grid_auto_load_max_variants

    def _set(self, override, value):
        override.node_grid_auto_load_max_variants = value
        override.save()

    def test_user_overrides_all(self):
        self._set(self.global_settings, 100)
        self._set(self.org_override, 200)
        self._set(self.lab_override, 300)
        self._set(self.user_override, 400)
        self.assertEqual(self._resolved(), 400)

    def test_lab_overrides_org_and_global(self):
        self._set(self.global_settings, 100)
        self._set(self.org_override, 200)
        self._set(self.lab_override, 300)
        self.assertEqual(self._resolved(), 300)

    def test_org_overrides_global(self):
        self._set(self.global_settings, 100)
        self._set(self.org_override, 200)
        self.assertEqual(self._resolved(), 200)

    def test_global_used_when_only_level_set(self):
        self._set(self.global_settings, 100)
        self.assertEqual(self._resolved(), 100)

    def test_falls_back_to_settings_constant_when_all_null(self):
        # Global is seeded with the code default by migration 0191 - clear it to exercise the
        # all-null path that the view's `... or settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS` covers.
        self._set(self.global_settings, None)
        self.assertIsNone(self._resolved())
        max_variants = self._resolved() or settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS
        self.assertEqual(max_variants, settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS)
