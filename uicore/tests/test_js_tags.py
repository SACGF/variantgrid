import json
from types import SimpleNamespace

from django.template import Context, RequestContext, Template
from django.test import RequestFactory, SimpleTestCase, TestCase

from uicore.templatetags.js_tags import jsonify_for_js


class JsonifyTest(SimpleTestCase):
    HOSTILE = 'a"b\\c\n</script><!-- & \u2028'

    def test_script_output_is_json_with_no_markup(self):
        for value in [self.HOSTILE, {"k": self.HOSTILE}, [self.HOSTILE]]:
            with self.subTest(value=value):
                js = jsonify_for_js(value)
                self.assertEqual(json.loads(js), value)
                for c in '<>&\n\u2028':
                    self.assertNotIn(c, js)

    def test_code_json_is_html_escaped(self):
        html = Template("{% load js_tags %}{% code_json data %}").render(Context({"data": {"k": "<img src=x>"}}))
        self.assertNotIn("<img", html)
        self.assertIn("&lt;img src=x&gt;", html)


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
