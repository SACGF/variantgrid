from unittest import TestCase

from genes.hgvs import PHGVS


class PhgvsTest(TestCase):

    def test_parts(self):
        parts = PHGVS.parse("p.Asp1054Ilefs*5", override_is_confirmed_to=False)
        assert not parts.transcript
        assert parts.aa_from == "Asp"
        assert "1054" == parts.codon == "1054"
        assert parts.aa_to == "Ile"
        assert parts.extra == "fs*5"
        assert parts.intron == False
        assert parts.is_confirmed == False
        assert parts.full_p_hgvs == "p.(Asp1054Ilefs*5)"

    def test_expand(self):
        parts = PHGVS.parse("p.Q1484*", override_is_confirmed_to=False)
        assert parts.full_p_hgvs == "p.(Gln1484*)"