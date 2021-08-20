from unittest import TestCase

from uicore.json.json_utils import JsonDiff, JsonDiffs
from uicore.json.validated_json import JsonMessages, ValidatedJson


class JsonUtilTests(TestCase):

    def test_validate_json(self):
        original = ValidatedJson({
            "normal": 3,
            "validated": ValidatedJson(
                [
                    1,
                    ValidatedJson(2, JsonMessages.info("This is the number 2"))
                ], JsonMessages.info("This is a list of numbers"))
        })
        serialized = original.serialize()
        deserialized = ValidatedJson.deserialize(serialized)
        self.assertEqual(original, deserialized)

        self.assertEqual(original.pure_json(), {"normal": 3, "validated": [1, 2]})


    def test_json_diff(self):
        diffs = JsonDiffs.differences(
            {"a": 3, "b": {"c": False}, 4: [1,2]},
            {"a": 4, "b": {"d": True}, 4: [1,2]}
        ).json_diffs
        print(diffs)
        self.assertEqual(len(diffs), 3)
        self.assertEqual(diffs[0].json_path_str, '["a"]')
        self.assertEqual(diffs[0].a, 3)
        self.assertEqual(diffs[0].b, 4)

        self.assertEqual(diffs[1].json_path_str, '["b"]["c"]')
        self.assertEqual(diffs[1].a, False)
        self.assertEqual(diffs[1].b, None)

        self.assertEqual(diffs[2].json_path_str, '["b"]["d"]')
        self.assertEqual(diffs[2].a, None)
        self.assertEqual(diffs[2].b, True)

        diffs = JsonDiffs.differences(
            ["A", "b"],
            ["a", "b", "c"]
        ).json_diffs
        print(diffs)