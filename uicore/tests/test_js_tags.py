import json

from django.test import SimpleTestCase

from uicore.templatetags.js_tags import jsonify


class JsonifyTest(SimpleTestCase):

    def test_string_is_a_js_string_literal(self):
        text = 'Say "hi"\\ \n</script>'
        js = jsonify(text)
        self.assertNotIn("</script>", js)
        self.assertEqual(json.loads(js), text)
