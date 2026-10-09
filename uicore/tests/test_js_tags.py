import json

from django.test import SimpleTestCase

from uicore.templatetags.js_tags import jsonify, jsonify_pretty

NASTY_TEXT = 'Say "hi"\\ \n`code` ${alert(1)} </SCRIPT><img src=x onerror=alert(1)> <!-- &amp;'


class JsonifyTest(SimpleTestCase):

    def _assert_safe_js_literal(self, js: str, value):
        for html_char in "<>&":
            self.assertNotIn(html_char, js)
        self.assertEqual(json.loads(js), value)

    def test_string_is_a_js_string_literal(self):
        self._assert_safe_js_literal(jsonify(NASTY_TEXT), NASTY_TEXT)

    def test_dict_is_a_js_literal(self):
        data = {"text": NASTY_TEXT, "</script>": [NASTY_TEXT]}
        self._assert_safe_js_literal(jsonify(data), data)

    def test_pretty_is_html_escaped(self):
        html = jsonify_pretty({"text": NASTY_TEXT})
        self.assertNotIn("<", html)
        self.assertIn("&lt;/SCRIPT&gt;", html)
