from django.test import TestCase
from classification.models.classification_patcher import patch_fuzzy_age


class ClassificationTestCaseAge(TestCase):

    def test_patch_fuzzy_age(self):
        self.assertEqual(patch_fuzzy_age("3"), "0-9")
        self.assertEqual(patch_fuzzy_age("17"), "10-19")
        self.assertEqual(patch_fuzzy_age("99"), "80+")
