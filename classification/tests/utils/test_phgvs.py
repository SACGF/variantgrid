from django.test import SimpleTestCase, TestCase
from hgvs.parsers import Parser

from genes.hgvs import PHGVS
from genes.hgvs.phgvs_canonical import canonical_p_hgvs, is_frameshift_in_disguise, short_p_hgvs


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


class CanonicalPHGVSTest(SimpleTestCase):
    """ ResolvedVariantInfo.p_hgvs form: source's accession kept, 3-letter, Ter, predicted brackets (#1383) """

    def test_canonical(self):
        for raw, expected in [
            ("p.G12D", "p.(Gly12Asp)"),
            ("NP_000483.3(CFTR):p.Gly480Cys", "NP_000483.3:p.(Gly480Cys)"),
            ("NP_000483.3:p.?", "NP_000483.3:p.?"),
            ("p.Met246Ilefs*14", "p.(Met246IlefsTer14)"),
            ("p.Trp26X", "p.(Trp26Ter)"),  # biocommons alone reads X as Xaa
            ("p.Arg97fsX5", "p.(Arg97fsTer5)"),
            ("ENSP00000371342.4:p.%3D", "ENSP00000371342.4:p.(=)"),
            ("p.=", "p.(=)"),
            ("p.?", "p.?"),
            ("p.(?)", "p.?"),
            ("p.Lys550_Lys558del", "p.(Lys550_Lys558del)"),
            ("c.35G>A p.(G12D)", "p.(Gly12Asp)"),
        ]:
            with self.subTest(raw=raw):
                self.assertEqual(expected, canonical_p_hgvs(raw))

    def test_not_p_hgvs(self):
        for raw in ["junk", "p.BLAGAMUDO", "NP_001139133.1:n.123+1G>A", "p.Glx2Lsy", "", None]:
            with self.subTest(raw=raw):
                self.assertIsNone(canonical_p_hgvs(raw))

    def test_short(self):
        self.assertEqual("p.(M246Ifs*14)", short_p_hgvs("NP_000483.3:p.(Met246IlefsTer14)"))

    def test_frameshift_in_disguise(self):
        parser = Parser()
        disguised = parser.parse("NP_1.1:p.Ser1321_Val2843delinsIleLeuTer")
        self.assertTrue(is_frameshift_in_disguise(disguised, protein_length=2843))
        self.assertFalse(is_frameshift_in_disguise(disguised, protein_length=3418))
        self.assertFalse(is_frameshift_in_disguise(parser.parse("NP_1.1:p.Cys28delinsTrpVal"), protein_length=100))
