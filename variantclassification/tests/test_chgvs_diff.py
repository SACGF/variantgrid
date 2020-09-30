from django.test.testcases import TestCase
from genes.hgvs import CHGVS, CHGVSDiff


def diff(c1: str, c2: str) -> CHGVSDiff:
    return CHGVS(c1).diff(CHGVS(c2))


# noinspection SpellCheckingInspection,SpellCheckingInspection
class CHGVSDiffTests(TestCase):

    def test_insert_gene(self):
        cd = diff('NM_000022.3:c.22G>A', 'NM_000022.3(ADA):c.22G>A')
        self.assertEqual(cd, CHGVSDiff.SAME)

    def test_same(self):
        cd = diff('NM_000071.2(CBS):c.919G>A', 'NM_000071.2(CBS):c.919G>A')
        self.assertEqual(cd, CHGVSDiff.SAME)

    def test_gene_change(self):
        cd = diff('NM_000071.2(CBS):c.919G>A', 'NM_000071.2(BOO):c.919G>A')
        self.assertEqual(cd, CHGVSDiff.DIFF_GENE)

    def test_transcript_change(self):
        cd = diff('NM_000099.2(CBS):c.919G>A', 'NM_000071.2(CBS):c.919G>A')
        self.assertEqual(cd, CHGVSDiff.DIFF_TRANSCRIPT_ID)

    def test_transcript_ver_change(self):
        cd = diff('NM_000071.2(CBS):c.919G>A', 'NM_000071.3(CBS):c.919G>A')
        self.assertEqual(cd, CHGVSDiff.DIFF_TRANSCRIPT_VER)

        # no version compared to a version should not raise a difference
        # otherwise every record from a client (that doesn't provide transcript versions)
        # will generate a warning
        cd = diff('NM_000071(CBS):c.919G>A', 'NM_000071.3(CBS):c.919G>A')
        self.assertEqual(cd, CHGVSDiff.SAME)

    def test_expanded(self):
        cd = diff('NM_001083961.1(WDR62):c.4205_4208del', 'NM_001083961.1(WDR62):c.4205_4208delTGCC')
        self.assertEqual(cd, CHGVSDiff.DIFF_RAW_CGVS_EXPANDED)

        cd = diff('NM_001083961.1(WDR62):c.4205_4208del4', 'NM_001083961.1(WDR62):c.4205_4208delTGCC')
        self.assertEqual(cd, CHGVSDiff.DIFF_RAW_CGVS_EXPANDED)

        # size has changed
        cd = diff('NM_001083961.1(WDR62):c.4205_4208del4', 'NM_001083961.1(WDR62):c.4205_4208delT')
        self.assertEqual(cd, CHGVSDiff.DIFF_RAW_CGVS)

        cd = diff('NM_001083961.1(WDR62):c.4205_4208del4', 'NM_001083961.1(WDR62):c.4205_4208del5')
        self.assertEqual(cd, CHGVSDiff.DIFF_RAW_CGVS)

        cd = diff('NM_001083961.1(WDR62):c.4205_4208delTGCC', 'NM_001083961.1(WDR62):c.4205_4208delTGGG')
        self.assertEqual(cd, CHGVSDiff.DIFF_RAW_CGVS)

        cd = diff('NM_001083961.1(WDR62):c.444205_4208del', 'NM_001083961.1(WDR62):c.4205_4208delTGCC')
        self.assertEqual(cd, CHGVSDiff.DIFF_RAW_CGVS)

    def test_multiple(self):
        cd = diff('NM_000081.2(CBS):c.919G>A', 'NM_000071.2(XXX):c.919G>C')
        self.assertEqual(cd, CHGVSDiff.DIFF_GENE | CHGVSDiff.DIFF_TRANSCRIPT_ID | CHGVSDiff.DIFF_RAW_CGVS)

    def test_multiple_2(self):
        # note we only complain about the gene if the gene was provided
        cd = diff('NM_006086.3:c.1012delC', 'NM_006086.3(TUBB3):c.1012delA')
        self.assertEqual(cd, CHGVSDiff.DIFF_RAW_CGVS)
