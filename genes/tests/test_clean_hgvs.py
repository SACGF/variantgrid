from django.test import TestCase

from genes.hgvs import HGVSMatcher


class TestCleanHGVS(TestCase):
    """ HGVSMatcher.clean_hgvs is now a thin wrapper over cdot.hgvs.clean_hgvs (a pure
        string operation - no genome build / converter / DB needed). """

    def test_clean_hgvs(self):
        test_cases = [
            ("M_206933.3(USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("NM_206933.3(USH2A)::c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("M_206933.3USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("m_206933.3USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("nm_206933.3(USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("m_206933.3(USH2A)::c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("m_206933.3(USH2A):c.4298G>A", "NM_206933.3(USH2A):c.4298G>A"),
            ("c_000007.13:g.117199563G>T", "NC_000007.13:g.117199563G>T"),
            ("c_000007.13::g.117199563G>T", "NC_000007.13:g.117199563G>T"),
            ("nc_000007.13:g.117199563G>T", "NC_000007.13:g.117199563G>T"),
            ("R_001566.1(TERC):n.427_428insC", "NR_001566.1(TERC):n.427_428insC"),
            ("r_001566.1(TERC)::n.427_428insC", "NR_001566.1(TERC):n.427_428insC"),
        ]

        for input_hgvs, expected_result in test_cases:
            result = HGVSMatcher.clean_hgvs(input_hgvs)
            self.assertEqual(result[0], expected_result, msg=f"Cleaning: {input_hgvs}")

    def test_clean_hgvs_corruption_canaries(self):
        """ Boundary canary for the cdot dependency: inputs where a cleaning
            regression would *silently corrupt* VG variant search. This is NOT a
            characterisation of cdot's cleaning (that's cdot's own suite) - it's the
            handful of VG-critical invariants that must survive any cdot version we
            pin, caught in VG CI rather than in production. """
        INVARIANTS = [
            # A gene symbol that contains a mutation-type substring (INS/DEL/DUP/INV)
            # must not be lower-cased into a different/invalid symbol. Old VG turned
            # INSR -> insR, breaking resolution for a real, common gene.
            ("NM_000208.4(INSR):c.215A>G", "NM_000208.4(INSR):c.215A>G"),
            # Balanced uncertain-range brackets are valid HGVS and must survive (old VG
            # stripped them, corrupting the variant).
            ("NM_004006.2(DMD):c.(4071+1_4072-1)_(5154+1_5155-1)del",
             "NM_004006.2(DMD):c.(4071+1_4072-1)_(5154+1_5155-1)del"),
            # Trailing protein annotation (common when pasting from reports) is removed
            # so the c.HGVS resolves.
            ("NM_000059.4:c.316+5G>A p.Arg106*", "NM_000059.4:c.316+5G>A"),
        ]
        for input_hgvs, expected in INVARIANTS:
            cleaned, _messages = HGVSMatcher.clean_hgvs(input_hgvs)
            self.assertEqual(cleaned, expected, msg=f"clean_hgvs({input_hgvs!r})")
