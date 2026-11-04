from django.test import SimpleTestCase

from genes.hgvs.hgvs import HGVSDisplay, abbreviate_hgvs_sequences
from genes.hgvs.phgvs_canonical import abbreviate_p_hgvs_sequences, short_p_hgvs
from snpdb.models.models_variant import abbreviate_variant_coordinate


class AbbreviateHGVSSequencesTest(SimpleTestCase):

    def test_abbreviate_hgvs_sequences(self):
        twenty = "A" * 20
        long_seq = "ACGT" * 5531
        self.assertEqual(abbreviate_hgvs_sequences(f"c.1_2ins{twenty}"), f"c.1_2ins{twenty}")
        self.assertEqual(abbreviate_hgvs_sequences(f"c.1_2ins{long_seq}"), "c.1_2ins<22,124 bp>")
        self.assertEqual(abbreviate_hgvs_sequences(f"c.5dup{twenty}C"), "c.5dup<21 bp>")
        # only the long side of a delins is shortened
        self.assertEqual(abbreviate_hgvs_sequences(f"c.1_30del{long_seq}insTG"), "c.1_30del<22,124 bp>insTG")
        self.assertEqual(abbreviate_hgvs_sequences(f"c.1_2delAGins{long_seq}"), "c.1_2delAGins<22,124 bp>")
        self.assertIsNone(abbreviate_hgvs_sequences(None))

    def test_hgvs_display(self):
        hgvs = HGVSDisplay.parse(dna_hgvs=f"NM_000059.4(BRCA2):c.1234_1235ins{'T' * 100}")
        self.assertTrue(hgvs.is_abbreviated)
        self.assertEqual(hgvs.dna_nomen_display, "c.1234_1235ins<100 bp>")
        self.assertFalse(HGVSDisplay.parse(dna_hgvs="NM_000059.4(BRCA2):c.1234G>A").is_abbreviated)


class AbbreviatePHGVSSequencesTest(SimpleTestCase):

    def test_abbreviate_p_hgvs_sequences(self):
        one_letter = "PTIAGTLAGKLRTKLKLYGSK"  # 21 residues
        self.assertEqual(abbreviate_p_hgvs_sequences(f"p.(Q1756delins{one_letter})"), "p.(Q1756delins<21 aa>)")
        self.assertEqual(abbreviate_p_hgvs_sequences(f"p.(Q1756delins{one_letter[:20]})"), f"p.(Q1756delins{one_letter[:20]})")
        three_letter = "ProThrIle" * 7
        self.assertEqual(abbreviate_p_hgvs_sequences(f"NP_000050.3:p.(Gln1756delins{three_letter})"),
                         "NP_000050.3:p.(Gln1756delins<21 aa>)")
        self.assertEqual(abbreviate_p_hgvs_sequences(f"p.(Q1756_A1757ins{'G' * 1000}*)"), "p.(Q1756_A1757ins<1,001 aa>)")
        self.assertEqual(abbreviate_p_hgvs_sequences("p.(M246Ifs*14)"), "p.(M246Ifs*14)")
        self.assertIsNone(abbreviate_p_hgvs_sequences(None))

    def test_short_p_hgvs_abbreviates(self):
        three_letter = "ProThrIle" * 7
        self.assertEqual(short_p_hgvs(f"NP_000050.3:p.(Gln1756delins{three_letter})"), "p.(Q1756delins<21 aa>)")


class AbbreviateVariantCoordinateTest(SimpleTestCase):

    def test_abbreviate_variant_coordinate(self):
        long_alt = "T" + "ACGT" * 5531
        self.assertEqual(abbreviate_variant_coordinate(f"17:43057062 T>{long_alt}"), "17:43057062 T>TAC...CGT")
        self.assertEqual(abbreviate_variant_coordinate(f"17:43057062 {long_alt}>T"), "17:43057062 TAC...CGT>T")
        self.assertEqual(abbreviate_variant_coordinate("17:43057062 T>TTTGAGCCCTTTGAGCCCTT"), "17:43057062 T>TTTGAGCCCTTTGAGCCCTT")
        self.assertEqual(abbreviate_variant_coordinate("17:43057061-43057062 <DUP>"), "17:43057061-43057062 <DUP>")
        self.assertIsNone(abbreviate_variant_coordinate(None))
