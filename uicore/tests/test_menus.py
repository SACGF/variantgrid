from collections import defaultdict
from unittest.mock import patch

from django.contrib import messages
from django.contrib.auth.models import User
from django.contrib.messages.storage import default_storage
from django.http import HttpResponse
from django.test import RequestFactory, SimpleTestCase, TestCase
from django.urls import get_resolver, reverse

from uicore.menus import current_menu
from uicore.page_frame import menu_html
from variantgrid.menus import MENUS


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


class PageFrameTest(TestCase):

    def setUp(self):
        menu_html.cache_clear()

    def tearDown(self):
        menu_html.cache_clear()

    def _frame(self, url_name, user=None) -> dict:
        if user:
            self.client.force_login(user)
        return self.client.get(reverse('page_frame'), {'url_name': url_name}).json()

    def test_highlight_and_admin_only(self):
        with _visible_except():
            user_frame = self._frame('variant_tags', User.objects.create_user('menu_user'))
            self.assertRegex(user_frame['menu_sub_html'], r'id="submenu-variant_tags"\s+class="nav-link active')
            self.assertRegex(user_frame['menu_main_html'], r'id="menu-top-variants"\s+class="nav-link active')
            self.assertNotIn('submenu-liftover_runs', user_frame['menu_sub_html'])

            superuser_frame = self._frame('variant_tags', User.objects.create_superuser('menu_admin'))
            self.assertIn('submenu-liftover_runs', superuser_frame['menu_sub_html'])

    def test_anonymous_gets_empty_frame(self):
        frame = self._frame('variant_tags')
        self.assertEqual(frame['menu_main_html'], '')
        self.assertEqual(frame['username'], '')

    def test_messages_from_previous_request_shown_once(self):
        user = User.objects.create_user('menu_messages')
        self.client.force_login(user)
        request = RequestFactory().get('/')
        # Store a message the way a view that redirects does, in this client's session
        session = self.client.session
        request.session = session
        storage = default_storage(request)
        storage.add(messages.INFO, "Saved the thing")
        response = HttpResponse()
        storage.update(response)
        session.save()
        self.client.cookies.update(response.cookies)

        self.assertIn("Saved the thing", self._frame('variant_tags')['messages_html'])
        self.assertEqual(self._frame('variant_tags')['messages_html'], '')
