import os
import unittest

from django.conf import settings
from django.test import TestCase

from ontology.management.commands import ontology_import


class Test(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    def testLoadData(self):
        small_owl = os.path.join(settings.BASE_DIR, "ontology", "tests", "test_data", "small.owl")
        ontology_import.load_hpo(small_owl, True)


if __name__ == "__main__":
    unittest.main()
