"""
Runs the shared HGVS corpus (genes/tests/test_data/hgvs_corpus.tsv) through HGVSComponents.

The corpus deliberately contains malformed input, because parsing is lenient by design -
it has to cope with whatever is already in the database. So these assert totality and
self-consistency rather than correctness of any individual parse.
"""
from django.test import TestCase

from genes.hgvs import HGVSComponents, HGVSDiff, HGVSDisplay
from genes.tests.utils.hgvs_corpus import all_hgvs, load_hgvs_corpus


class HGVSCorpusTests(TestCase):

    def test_corpus_loads(self):
        by_category = load_hgvs_corpus()
        self.assertGreater(len(all_hgvs()), 400)
        for category in ("refseq_c", "ensembl_c", "with_gene_symbol", "malformed", "lrg"):
            self.assertIn(category, by_category)

    def test_full_hgvs_round_trips(self):
        """ Parsing is total - every input comes back out whole, however malformed """
        for hgvs in all_hgvs():
            with self.subTest(hgvs=hgvs):
                self.assertEqual(str(HGVSComponents(hgvs)), hgvs)

    def test_rewriters_and_json_never_raise(self):
        for hgvs in all_hgvs():
            with self.subTest(hgvs=hgvs):
                components = HGVSComponents(hgvs)
                self.assertIsNotNone(components.sort_str)
                self.assertIsNotNone(components.without_gene_symbol_str)
                self.assertIsNotNone(components.without_transcript_version)
                self.assertIsNotNone(components.with_transcript_version(9))
                self.assertIsNotNone(components.with_gene_symbol("BRCA1"))
                self.assertIn("full", HGVSDisplay(components).to_json())

    def test_sorting_is_a_total_order(self):
        """ Order must depend on the values alone, never on the order they arrived in -
            a key where one entry is a prefix of another is easy to get wrong """
        corpus = all_hgvs()
        forwards = [str(d) for d in sorted(HGVSDisplay.parse(hgvs) for hgvs in corpus)]
        backwards = [str(d) for d in sorted(HGVSDisplay.parse(hgvs) for hgvs in reversed(corpus))]
        self.assertEqual(forwards, backwards)

    def test_diff_against_self_is_same(self):
        for hgvs in all_hgvs():
            with self.subTest(hgvs=hgvs):
                self.assertEqual(HGVSComponents(hgvs).diff(HGVSComponents(hgvs)), HGVSDiff.SAME)

    def test_diff_is_symmetric(self):
        corpus = all_hgvs()
        for a, b in zip(corpus, corpus[1:]):
            with self.subTest(a=a, b=b):
                self.assertEqual(HGVSComponents(a).diff(HGVSComponents(b)),
                                 HGVSComponents(b).diff(HGVSComponents(a)))

    def test_equal_components_hash_equal(self):
        for hgvs in all_hgvs():
            with self.subTest(hgvs=hgvs):
                a, b = HGVSComponents(hgvs), HGVSComponents(hgvs)
                self.assertEqual(a, b)
                self.assertEqual(len({a, b}), 1)
