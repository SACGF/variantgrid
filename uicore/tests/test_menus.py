from collections import defaultdict
from unittest.mock import patch

from django.contrib.auth.models import AnonymousUser, User
from django.template import RequestContext, Template
from django.test import RequestFactory, SimpleTestCase, TestCase
from django.urls import get_resolver, resolve, reverse

from uicore.menus import MENUS, current_menu


def _visible_except(*hidden):
    visible = defaultdict(lambda: True)
    for url_name in hidden:
        visible[url_name] = False
    return patch('uicore.menus.get_visible_url_names', return_value=visible)


class MenuRegistryTest(SimpleTestCase):

    def test_url_names_resolve(self):
        url_names = {k for k in get_resolver().reverse_dict.keys() if isinstance(k, str)}
        for menu in MENUS:
            declared = list(menu.pages) + ([menu.url_name] if menu.url_name else [])
            for item in menu.items:
                declared.extend(item.pages)
                if not item.href:
                    declared.append(item.url_name)
            unknown = set(declared) - url_names
            self.assertFalse(unknown, f"menu '{menu.key}' names unknown url names")

    def test_detail_page_menu(self):
        with _visible_except():
            self.assertEqual(current_menu('view_allele').key, 'variants')
            self.assertEqual(current_menu('view_liftover_run').key, 'variants')
            self.assertIsNone(current_menu('no_such_page'))

    def test_liftover_moves_to_settings_without_variants_menu(self):
        with _visible_except('variant_tags', 'variants'):
            self.assertEqual(current_menu('liftover_runs').key, 'settings')
            self.assertEqual(current_menu('view_liftover_run').key, 'settings')

    def test_page_under_hidden_item_keeps_menu(self):
        """ A superuser can reach a page whose item this deployment hides - it still gets its menu """
        with _visible_except('clinvar_key_summary'):
            self.assertEqual(current_menu('clinvar_export').key, 'classifications')


class MenuBarSubTest(TestCase):

    def _render(self, url_name, user) -> str:
        request = RequestFactory().get(reverse(url_name))
        request.resolver_match = resolve(request.path)
        request.user = user
        return Template("{% load ui_menus %}{% menu_bar_sub %}").render(RequestContext(request))

    def test_highlight_and_admin_only(self):
        with _visible_except():
            html = self._render('variant_tags', AnonymousUser())
            self.assertRegex(html, r'id="submenu-variant_tags"\s+class="nav-link active')
            self.assertNotIn('submenu-liftover_runs', html)

            superuser = User.objects.create_superuser('menu_admin')
            self.assertIn('submenu-liftover_runs', self._render('variant_tags', superuser))
