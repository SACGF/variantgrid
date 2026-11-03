from django.test import SimpleTestCase

from genes.hgvs.hgvs import HGVSDisplay, abbreviate_hgvs_sequences


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
