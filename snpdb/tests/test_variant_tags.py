"""
snpdb.templatetags.variant_tags: the friendly size a structural variant shows on the variant page.
"""
from typing import Optional

from django.test import SimpleTestCase

from snpdb.models.models_variant import Sequence, Variant
from snpdb.templatetags.variant_tags import sv_size


def _variant(alt: str, svlen: Optional[int]) -> Variant:
    return Variant(alt=Sequence(seq=alt), svlen=svlen)


class SVSizeTest(SimpleTestCase):
    def test_kind_and_size(self):
        cases = [
            ("<DEL>", -11_000_000, "11 Mb deletion"),
            ("<DUP>", 1_500, "1.5 kb duplication"),
            ("<INV>", 55_123, "55.12 kb inversion"),
            ("<INS>", 50, "50 bp insertion"),
            ("<CNV>", 2_000_000, "2 Mb CNV"),
        ]
        for alt, svlen, expected in cases:
            with self.subTest(alt=alt):
                self.assertEqual(sv_size(_variant(alt, svlen)), expected)

    def test_no_size_without_a_symbolic_sv_alt(self):
        self.assertEqual(sv_size(_variant("A", None)), "")
        # Gene-level events are symbolic but have no length
        self.assertEqual(sv_size(_variant("<FUSION:HGNC:1100>", 0)), "")
