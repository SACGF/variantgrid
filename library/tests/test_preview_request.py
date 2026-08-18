from datetime import datetime

from django.test import TestCase

from genes.models_enums import AnnotationConsortium
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
    """ PreviewData holds unhashable collection fields, so __hash__ is hand written. """

    def test_equal_objects_have_same_hash(self):
        pd1 = PreviewData(category="Gene", identifier="BRCA1", title="BRCA1 gene")
        pd2 = PreviewData(category="Gene", identifier="BRCA1", title="BRCA1 gene")
        self.assertEqual(hash(pd1), hash(pd2))

    def test_different_objects_have_different_hash(self):
        pd1 = PreviewData(category="Gene", identifier="BRCA1")
        pd2 = PreviewData(category="Gene", identifier="BRCA2")
        self.assertNotEqual(hash(pd1), hash(pd2))

    def test_summary_extra_list_hashed_by_key_and_value(self):
        def _pd(value):
            return PreviewData(category="Gene", identifier="BRCA1",
                               summary_extra=[PreviewKeyValue(key="count", value=value)])

        self.assertEqual(hash(_pd(42)), hash(_pd(42)))
        self.assertNotEqual(hash(_pd(42)), hash(_pd(43)))

    def test_annotation_consortia_set_hashed_regardless_of_order(self):
        def _pd(*consortia):
            return PreviewData(category="Gene", identifier="BRCA1",
                               annotation_consortia=set(consortia))

        both_orders = (
            _pd(AnnotationConsortium.REFSEQ, AnnotationConsortium.ENSEMBL),
            _pd(AnnotationConsortium.ENSEMBL, AnnotationConsortium.REFSEQ),
        )
        self.assertEqual(hash(both_orders[0]), hash(both_orders[1]))
        self.assertNotEqual(hash(both_orders[0]), hash(_pd(AnnotationConsortium.ENSEMBL)))
