from django.test import SimpleTestCase

from genes.hgvs.phgvs_canonical import canonical_p_hgvs, short_p_hgvs


class PHGVSCanonicalTest(SimpleTestCase):

    def test_initiation_codon(self):
        """ A change to Met1 stops translation starting - HGVS says never write it as p.Met1Thr (#1383) """
        expected = {
            "NP_055635.3:p.Met1Thr": "NP_055635.3:p.(Met1?)",  # ClinGen
            "p.M1T": "p.(Met1?)",
            "p.Met1_Ala4del": "p.(Met1?)",
            "p.Met1ValfsTer5": "p.(Met1?)",
            "ENSP1.1:p.Met1?": "ENSP1.1:p.(Met1?)",  # VEP
            "p.Met1ext-5": "p.(Met1ext-5)",  # a new upstream start
            "p.Met1_Thr6dup": "p.(Met1_Thr6dup)",
        }
        for raw, canonical in expected.items():
            with self.subTest(raw=raw):
                self.assertEqual(canonical, canonical_p_hgvs(raw))
        self.assertEqual("p.(M1?)", short_p_hgvs("NP_055635.3:p.(Met1?)"))
