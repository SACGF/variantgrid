from unittest import TestCase

from genes.hgvs import PHGVS


class PhgvsTest(TestCase):

    def test_parts(self):
        parts = PHGVS.parse("p.Asp1054Ilefs*5", override_is_confirmed_to=False)
        self.assertFalse(bool(parts.transcript))
        self.assertEqual(parts.aa_from, "Asp")
        self.assertEqual(parts.codon, "1054")
        self.assertEqual(parts.aa_to, "Ile")
        self.assertEqual(parts.extra, "fs*5")
        self.assertEqual(parts.intron, False)
        self.assertFalse(parts.is_confirmed, False)
        self.assertEqual(parts.full_p_hgvs, "p.(Asp1054Ilefs*5)")

    def test_expand(self):
        parts = PHGVS.parse("p.Q1484*", override_is_confirmed_to=False)
        self.assertEqual(parts.full_p_hgvs, "p.(Gln1484*)")

    def test_equal(self):
        p1 = PHGVS.parse("p.Pro81Leu", override_is_confirmed_to=False).without_transcript
        p2 = PHGVS.parse("NP_002993.1:p.(Pro81Leu)", override_is_confirmed_to=False).without_transcript
        self.assertEqual(p1, p2)
        self.assertEqual(hash(p1), hash(p2))

    def test_intron(self):
        p1 = PHGVS.parse("p.?")
        self.assertTrue(p1.intron)
        self.assertEqual(p1.full_p_hgvs, "p.?")

    def test_transcript(self):
        p1 = PHGVS.parse("NP_001269862.1:p.A808=", override_is_confirmed_to=False)
        self.assertEqual("NP_001269862.1", p1.transcript,)
        self.assertEqual("p.(Ala808=)", str(p1.without_transcript))

        # can't parse the p. but should at least be able to separate out the transcript
        p1 = PHGVS.parse("NP_001269862.1:p.BLAGAMUDO", override_is_confirmed_to=False)
        self.assertEqual("p.BLAGAMUDO", str(p1.without_transcript))
