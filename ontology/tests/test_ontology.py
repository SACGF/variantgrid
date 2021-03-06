import unittest

from django.test import TestCase

from ontology.tests.test_data_ontology import create_ontology_test_data


class Test(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    def testLoadData(self):
        create_ontology_test_data()


if __name__ == "__main__":
    unittest.main()
