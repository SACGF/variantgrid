from datetime import datetime

from django.test import TestCase

from library.preview_request import PreviewData, PreviewKeyValue


class TestPreviewKeyValueValueStr(TestCase):
    def test_datetime_value_str_formats_instance(self):
        """value_str should format the datetime instance, not the class"""
        dt = datetime(2024, 6, 15)
        pkv = PreviewKeyValue(key="date", value=dt)
        self.assertEqual(pkv.value_str, "2024-06-15")

    def test_string_value_str(self):
        pkv = PreviewKeyValue(key="name", value="hello")
        self.assertEqual(pkv.value_str, "hello")

    def test_int_value_str(self):
        pkv = PreviewKeyValue(key="count", value=42)
        self.assertEqual(pkv.value_str, "42")


class TestPreviewDataHash(TestCase):
    def test_hash_is_stable(self):
        """__hash__ should return the same value on repeated calls"""
        pd = PreviewData(category="Gene", identifier="BRCA1")
        self.assertEqual(hash(pd), hash(pd))

    def test_equal_objects_have_same_hash(self):
        pd1 = PreviewData(category="Gene", identifier="BRCA1", title="BRCA1 gene")
        pd2 = PreviewData(category="Gene", identifier="BRCA1", title="BRCA1 gene")
        self.assertEqual(hash(pd1), hash(pd2))

    def test_different_objects_have_different_hash(self):
        pd1 = PreviewData(category="Gene", identifier="BRCA1")
        pd2 = PreviewData(category="Gene", identifier="BRCA2")
        self.assertNotEqual(hash(pd1), hash(pd2))

    def test_hash_usable_in_set(self):
        """PreviewData instances should be usable in a set (requires stable hashing)"""
        pd = PreviewData(category="Gene", identifier="BRCA1")
        s = {pd}
        self.assertIn(pd, s)
