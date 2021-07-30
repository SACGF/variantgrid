from unittest import TestCase
import json
from classification.json_utils import ValidatedJson, JsonMessages


class JsonUtilTests(TestCase):

    def test_validate_json(self):
        original = ValidatedJson({
            "normal": 3,
            "validated": ValidatedJson(
                 [
                     1,
                     ValidatedJson(2, JsonMessages.info("This is the number 2"))
                 ]
            , JsonMessages.info("This is a list of numbers"))
        })
        serialized = original.serialize()
        print(json.dumps(serialized))
        deserialized = ValidatedJson.deserialize(serialized)
        self.assertEqual(original, deserialized)
