import json

from django.test import SimpleTestCase

from uicore.templatetags.js_tags import json_data


class JsonDataTest(SimpleTestCase):

    def test_non_finite_as_null(self):
        html = json_data("x-data", a=float("nan"), b=[1.5, float("inf")])
        json_str = html.removeprefix('<script id="x-data" type="application/json">').removesuffix("</script>")
        self.assertEqual(json.loads(json_str), {"a": None, "b": [1.5, None]})
