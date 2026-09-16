from collections import defaultdict

from django.apps import apps
from django.conf import settings
from django.contrib.auth.models import User
from django.template import Context, Template
from django.test import TestCase, override_settings
from django.urls.resolvers import get_resolver

from snpdb.models.models_user_settings import UserSettingsOverride
from variantgrid.perm_path import get_visible_url_names
from variantgrid.tips import _load_tips, _setting_enabled, get_tips


class TipsDataTestCase(TestCase):
    """ Keeps tips.csv honest as URLs get renamed and settings come and go """

    def test_urls_are_registered_names(self):
        registered = {k for k in get_resolver(None).reverse_dict.keys() if isinstance(k, str)}
        unknown = {tip.url for tip in _load_tips()} - registered
        self.assertFalse(unknown, f"tips.csv URL names not registered: {sorted(unknown)}")

    def test_apps_are_installed(self):
        installed = {app.label for app in apps.get_app_configs()}
        unknown = {tip.app for tip in _load_tips()} - installed
        self.assertFalse(unknown, f"tips.csv apps not installed: {sorted(unknown)}")

    def test_settings_enabled_paths_resolve(self):
        for path in sorted({tip.settings_enabled for tip in _load_tips()} - {""}):
            with self.subTest(path=path):
                value = settings
                for step in path.split("."):
                    if value is None:
                        break  # Feature switched off entirely (eg CLINVAR_EXPORT = None)
                    if isinstance(value, dict):
                        self.assertIn(step, value, f"'{path}' - no key '{step}'")
                        value = value[step]
                    else:
                        self.assertTrue(hasattr(value, step), f"'{path}' - no setting '{step}'")
                        value = getattr(value, step)


class SettingEnabledTestCase(TestCase):
    def test_blank_path_is_enabled(self):
        self.assertTrue(_setting_enabled(""))

    @override_settings(SOMALIER={"enabled": True})
    def test_dict_path(self):
        self.assertTrue(_setting_enabled("SOMALIER.enabled"))

    @override_settings(SOMALIER={"enabled": False})
    def test_dict_path_disabled(self):
        self.assertFalse(_setting_enabled("SOMALIER.enabled"))

    @override_settings(CLINVAR_EXPORT=None)
    def test_path_through_none(self):
        self.assertFalse(_setting_enabled("CLINVAR_EXPORT.mode"))

    def test_unknown_setting_is_disabled(self):
        self.assertFalse(_setting_enabled("NO_SUCH_SETTING"))


@override_settings(TIPS_ENABLED=True)
class GetTipsTestCase(TestCase):
    """ get_tips is cached, so anything that changes settings clears it first """

    def setUp(self):
        get_tips.cache_clear()
        self.addCleanup(get_tips.cache_clear)

    def test_tips_enabled_off_returns_nothing(self):
        with override_settings(TIPS_ENABLED=False):
            get_tips.cache_clear()
            self.assertEqual([], get_tips())

    def test_hidden_urls_are_excluded(self):
        self.assertIn("analysis", {tip.url for tip in get_tips()})

        url_name_register = defaultdict(lambda: True, settings.URLS_NAME_REGISTER)
        url_name_register["analysis"] = False
        with override_settings(URLS_NAME_REGISTER=url_name_register):
            get_tips.cache_clear()
            get_visible_url_names.cache_clear()
            self.addCleanup(get_visible_url_names.cache_clear)
            self.assertNotIn("analysis", {tip.url for tip in get_tips()})

    def test_settings_enabled_gates_tips(self):
        with override_settings(MME_ENABLED=False):
            get_tips.cache_clear()
            self.assertEqual([], [t for t in get_tips() if t.settings_enabled == "MME_ENABLED"])

        with override_settings(MME_ENABLED=True):
            get_tips.cache_clear()
            self.assertTrue([t for t in get_tips() if t.settings_enabled == "MME_ENABLED"])

    def test_somalier_gates_tips(self):
        with override_settings(SOMALIER={**settings.SOMALIER, "enabled": False}):
            get_tips.cache_clear()
            self.assertEqual([], [t for t in get_tips() if t.settings_enabled == "SOMALIER.enabled"])

        with override_settings(SOMALIER={**settings.SOMALIER, "enabled": True}):
            get_tips.cache_clear()
            self.assertTrue([t for t in get_tips() if t.settings_enabled == "SOMALIER.enabled"])

    def test_app_filter(self):
        analysis_tips = get_tips(app="analysis")
        self.assertTrue(analysis_tips)
        self.assertEqual({"analysis"}, {tip.app for tip in analysis_tips})
        self.assertLess(len(analysis_tips), len(get_tips()))


@override_settings(TIPS_ENABLED=True)
class TipBoxTagTestCase(TestCase):
    def setUp(self):
        get_tips.cache_clear()
        self.addCleanup(get_tips.cache_clear)
        self.user = User.objects.create_user("tips_test_user")
        self.template = Template("{% load tips_tags %}{% tip_box %}")

    def _render(self):
        return self.template.render(Context({"user": self.user})).strip()

    def test_renders_tip_for_user(self):
        self.assertIn("tip-box", self._render())

    def test_renders_nothing_when_user_turned_tips_off(self):
        UserSettingsOverride.objects.update_or_create(user=self.user, defaults={"show_tips": False})
        self.assertEqual("", self._render())

    def test_renders_nothing_when_tips_disabled(self):
        with override_settings(TIPS_ENABLED=False):
            get_tips.cache_clear()
            self.assertEqual("", self._render())
