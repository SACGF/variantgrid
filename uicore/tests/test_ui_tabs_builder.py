from types import SimpleNamespace

from django.template import RequestContext, Template
from django.test import RequestFactory, TestCase


class TabEmbeddedAdminOnlyTest(TestCase):
    TEMPLATE = ('{% load ui_tabs_builder %}'
                '{% ui_register_tab_embedded tab_set="t" label="Admin" admin_only=True %}secret{% end_ui_register_tab_embedded %}'
                '{% ui_render_tabs tab_set="t" %}')

    def _render(self, is_superuser: bool) -> str:
        request = RequestFactory().get("/")
        request.user = SimpleNamespace(is_superuser=is_superuser)
        return Template(self.TEMPLATE).render(RequestContext(request))

    def test_hidden_from_non_superuser(self):
        self.assertNotIn("secret", self._render(is_superuser=False))

    def test_shown_to_superuser(self):
        self.assertIn("secret", self._render(is_superuser=True))
