from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, TrioNode
from analysis.models.enums import TrioInheritance
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild, Variant, BuiltInFilters
from snpdb.models.models_cohort import CohortGenotype, CohortGenotypeCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestTrioNodeInheritance(TestCase):
    """
    End-to-end tests for TrioNode inheritance filters.

    Creates real Variant + CohortGenotype records and verifies that the correct
    variants are returned (and wrong ones are excluded) for each inheritance mode.

    Trio sample layout (from create_fake_trio → create_fake_cohort):
        packed_field_index 0 → proband
        packed_field_index 1 → mother  (mother_affected=True)
        packed_field_index 2 → father  (father_affected=False)

    samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='testuser_trio')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        # mother_affected=True, father_affected=False
        cls.trio = create_fake_trio(user, cls.grch37)
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cls.trio.cohort)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        # One variant per inheritance pattern.  samples_zygosity index: [proband, mother, father]
        #
        # Recessive: proband=HOM_ALT, both parents=HET
        cls.recessive_v = slowly_create_test_variant("3", 1000, "A", "T", cls.grch37)
        cls._make_cg(cls.recessive_v, "OEE")

        # Denovo: proband=HET, both parents=HOM_REF
        cls.denovo_v = slowly_create_test_variant("3", 2000, "A", "T", cls.grch37)
        cls._make_cg(cls.denovo_v, "ERR")

        # Dominant (mother affected): proband=HET, mother=HET, father=HOM_REF
        cls.dominant_v = slowly_create_test_variant("3", 3000, "A", "T", cls.grch37)
        cls._make_cg(cls.dominant_v, "EER")

        # X-linked recessive: proband=HOM_ALT, mother=HET, father=any — on chrX
        cls.xlinked_v = slowly_create_test_variant("X", 1000, "A", "T", cls.grch37)
        cls._make_cg(cls.xlinked_v, "OER")

        # Control: all HOM_REF — should not match any inheritance mode
        cls.control_v = slowly_create_test_variant("3", 4000, "A", "T", cls.grch37)
        cls._make_cg(cls.control_v, "RRR")

        # Unknown-father variant: proband=HOM_ALT, mother=HET, father=UNKNOWN
        # Used to test require_zygosity behaviour.
        cls.unknown_father_v = slowly_create_test_variant("3", 5000, "A", "T", cls.grch37)
        cls._make_cg(cls.unknown_father_v, "OEU")

        # X-linked with unknown mother: proband=HOM_ALT, mother=UNKNOWN, father=any
        # Used to test XLR-branch require_zygosity in ALL_RECESSIVE.
        cls.xlinked_unknown_mother_v = slowly_create_test_variant("X", 2000, "A", "T", cls.grch37)
        cls._make_cg(cls.xlinked_unknown_mother_v, "OUR")

        # Proband-only variant: proband=HET, mother=HOM_REF, father=HOM_REF
        cls.proband_only_v = slowly_create_test_variant("3", 6000, "A", "T", cls.grch37)
        cls._make_cg(cls.proband_only_v, "ERR")

        # Mother-only variant: proband=HOM_REF, mother=HET, father=HOM_REF
        cls.mother_only_v = slowly_create_test_variant("3", 6100, "A", "T", cls.grch37)
        cls._make_cg(cls.mother_only_v, "RER")

        # Father-only variant: proband=HOM_REF, mother=HOM_REF, father=HET
        cls.father_only_v = slowly_create_test_variant("3", 6200, "A", "T", cls.grch37)
        cls._make_cg(cls.father_only_v, "RRE")

    @classmethod
    def _make_cg(cls, variant, samples_zygosity):
        n = len(samples_zygosity)
        CohortGenotype.objects.create(
            collection=cls.cgc,
            variant=variant,
            ref_count=samples_zygosity.count('R'),
            het_count=samples_zygosity.count('E'),
            hom_count=samples_zygosity.count('O'),
            samples_zygosity=samples_zygosity,
            samples_allele_depth=[20] * n,
            samples_allele_frequency=[100] * n,
            samples_read_depth=[30] * n,
            samples_genotype_quality=[30] * n,
            samples_phred_likelihood=[0] * n,
        )

    def _make_node(self, inheritance, **kwargs):
        return TrioNode.objects.create(
            analysis=self.analysis, trio=self.trio,
            inheritance=inheritance, **kwargs
        )

    def _filter_variants(self, node):
        """Apply node's inheritance filter against real Variant data; return matching PKs."""
        arg_q_dict = node._get_node_arg_q_dict()
        cgc = node.trio.cohort.cohort_genotype_collection
        qs = Variant.objects.annotate(**cgc.get_annotation_kwargs())
        alias = cgc.cohortgenotype_alias
        if alias in arg_q_dict:
            for q in arg_q_dict[alias].values():
                qs = qs.filter(q)
        if None in arg_q_dict:
            for q in arg_q_dict[None].values():
                qs = qs.filter(q)
        return set(qs.values_list('pk', flat=True))

    # ── Recessive ─────────────────────────────────────────────────────────────

    def test_recessive_matches_recessive_variant(self):
        node = self._make_node(TrioInheritance.RECESSIVE)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_recessive_excludes_other_variants(self):
        node = self._make_node(TrioInheritance.RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.dominant_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── Denovo ────────────────────────────────────────────────────────────────

    def test_denovo_matches_denovo_variant(self):
        node = self._make_node(TrioInheritance.DENOVO)
        self.assertIn(self.denovo_v.pk, self._filter_variants(node))

    def test_denovo_excludes_other_variants(self):
        node = self._make_node(TrioInheritance.DENOVO)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.dominant_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── Dominant ──────────────────────────────────────────────────────────────

    def test_dominant_matches_dominant_variant(self):
        # trio has mother_affected=True, so dominant expects mother=HAS_VARIANT
        node = self._make_node(TrioInheritance.DOMINANT)
        self.assertIn(self.dominant_v.pk, self._filter_variants(node))

    def test_dominant_excludes_other_variants(self):
        node = self._make_node(TrioInheritance.DOMINANT)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── X-linked recessive ────────────────────────────────────────────────────

    def test_xlinked_matches_x_chromosome_variant(self):
        node = self._make_node(TrioInheritance.XLINKED_RECESSIVE)
        self.assertIn(self.xlinked_v.pk, self._filter_variants(node))

    def test_xlinked_excludes_autosomal_variants(self):
        node = self._make_node(TrioInheritance.XLINKED_RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.dominant_v.pk, ids)

    # ── require_zygosity ──────────────────────────────────────────────────────

    def test_require_zygosity_true_excludes_unknown_parental_zygosity(self):
        node = self._make_node(TrioInheritance.RECESSIVE, require_zygosity=True)
        self.assertNotIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_require_zygosity_false_allows_unknown_parental_zygosity(self):
        node = self._make_node(TrioInheritance.RECESSIVE, require_zygosity=False)
        self.assertIn(self.unknown_father_v.pk, self._filter_variants(node))

    # ── Validation ────────────────────────────────────────────────────────────

    def test_dominant_with_affected_mother_no_errors(self):
        # trio has mother_affected=True
        errors = TrioNode.get_trio_inheritance_errors(self.trio, TrioInheritance.DOMINANT)
        self.assertEqual(errors, [])

    def test_dominant_no_affected_parent_raises_error(self):
        self.trio.mother_affected = False
        self.trio.father_affected = False
        errors = TrioNode.get_trio_inheritance_errors(self.trio, TrioInheritance.DOMINANT)
        self.assertGreater(len(errors), 0)
        self.assertIn("affected parent", errors[0].lower())
        self.trio.mother_affected = True  # restore
        self.trio.save()

    def test_xlinked_affected_mother_raises_error(self):
        # trio has mother_affected=True — xlinked requires unaffected mother
        errors = TrioNode.get_trio_inheritance_errors(self.trio, TrioInheritance.XLINKED_RECESSIVE)
        self.assertGreater(len(errors), 0)
        self.assertIn("mother", errors[0].lower())

    # ── CompHet node type ─────────────────────────────────────────────────────

    def test_compound_het_requires_one_parent_input(self):
        node = self._make_node(TrioInheritance.COMPOUND_HET)
        self.assertEqual(node.min_inputs, 1)
        self.assertEqual(node.max_inputs, 1)

    def test_simple_modes_are_source_nodes(self):
        for mode in [TrioInheritance.RECESSIVE, TrioInheritance.DENOVO,
                     TrioInheritance.DOMINANT, TrioInheritance.XLINKED_RECESSIVE]:
            node = self._make_node(mode)
            self.assertEqual(node.max_inputs, 0, f"{mode} should be a source node")

    # ── Clone ─────────────────────────────────────────────────────────────────

    def test_clone_matches_same_variants(self):
        node = self._make_node(TrioInheritance.RECESSIVE)
        clone = node.save_clone()
        self.assertEqual(self._filter_variants(node), self._filter_variants(clone))

    # ── ALL_RECESSIVE (AR ∪ XLR) ──────────────────────────────────────────────

    def test_all_recessive_matches_autosomal_recessive_variant(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_all_recessive_matches_xlinked_variant(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE)
        self.assertIn(self.xlinked_v.pk, self._filter_variants(node))

    def test_all_recessive_excludes_unrelated_variants(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    def test_all_recessive_require_zygosity_excludes_unknown_father_on_autosome(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE, require_zygosity=True)
        self.assertNotIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_all_recessive_no_require_zygosity_includes_unknown_father_on_autosome(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE, require_zygosity=False)
        self.assertIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_all_recessive_require_zygosity_excludes_unknown_mother_on_xlr(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE, require_zygosity=True)
        self.assertNotIn(self.xlinked_unknown_mother_v.pk, self._filter_variants(node))

    def test_all_recessive_no_require_zygosity_includes_unknown_mother_on_xlr(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE, require_zygosity=False)
        self.assertIn(self.xlinked_unknown_mother_v.pk, self._filter_variants(node))

    def test_all_recessive_is_source_node(self):
        node = self._make_node(TrioInheritance.ALL_RECESSIVE)
        self.assertEqual(node.max_inputs, 0)

    def test_zygosity_table_all_recessive_two_line_cells(self):
        data = TrioNode.get_zygosity_table_data()
        self.assertIn('AR:', data[TrioInheritance.ALL_RECESSIVE]['mother'])
        self.assertIn('XLR:', data[TrioInheritance.ALL_RECESSIVE]['mother'])

    def test_zygosity_table_all_recessive_other_filters_mentions_chr_x(self):
        data = TrioNode.get_zygosity_table_data()
        self.assertIn('Chr X', data[TrioInheritance.ALL_RECESSIVE]['other_filters_mother'])

    # ── ANY_AFFECTED ──────────────────────────────────────────────────────────

    def test_any_affected_includes_variant_in_proband_only(self):
        # Trio fixture has mother_affected=True; proband is always affected.
        node = self._make_node(TrioInheritance.ANY_AFFECTED)
        self.assertIn(self.proband_only_v.pk, self._filter_variants(node))

    def test_any_affected_includes_variant_in_affected_mother(self):
        node = self._make_node(TrioInheritance.ANY_AFFECTED)
        self.assertIn(self.mother_only_v.pk, self._filter_variants(node))

    def test_any_affected_excludes_variant_in_unaffected_father_only(self):
        node = self._make_node(TrioInheritance.ANY_AFFECTED)
        self.assertNotIn(self.father_only_v.pk, self._filter_variants(node))

    def test_any_affected_collapses_when_no_affected_parent(self):
        # Flip the trio so only the proband is affected; variant in mother
        # alone should be excluded.
        self.trio.mother_affected = False
        self.trio.save()
        try:
            node = self._make_node(TrioInheritance.ANY_AFFECTED)
            ids = self._filter_variants(node)
            self.assertIn(self.proband_only_v.pk, ids)
            self.assertNotIn(self.mother_only_v.pk, ids)
        finally:
            self.trio.mother_affected = True
            self.trio.save()

    def test_any_affected_is_source_node(self):
        node = self._make_node(TrioInheritance.ANY_AFFECTED)
        self.assertEqual(node.max_inputs, 0)

    def test_any_affected_always_valid_no_errors(self):
        errors = TrioNode.get_trio_inheritance_errors(self.trio, TrioInheritance.ANY_AFFECTED)
        self.assertEqual(errors, [])

    def test_zygosity_table_any_affected_proband_has_variant(self):
        data = TrioNode.get_zygosity_table_data()
        proband_value = data[TrioInheritance.ANY_AFFECTED]['proband']
        self.assertTrue(proband_value)

    # ── Zygosity table "Other Filters" column ─────────────────────────────────

    def test_zygosity_table_xlinked_has_chr_x_other_filter(self):
        data = TrioNode.get_zygosity_table_data()
        self.assertEqual(data[TrioInheritance.XLINKED_RECESSIVE]['other_filters_mother'], "Chr X only")

    def test_zygosity_table_recessive_has_no_other_filter_keys(self):
        data = TrioNode.get_zygosity_table_data()
        keys = [k for k in data[TrioInheritance.RECESSIVE] if k.startswith('other_filters_')]
        self.assertEqual(keys, [])

    def test_zygosity_table_compound_het_has_gene_constraint_in_other_filters(self):
        data = TrioNode.get_zygosity_table_data()
        self.assertIn("gene", data[TrioInheritance.COMPOUND_HET]['other_filters_mother'])

    def test_zygosity_table_compound_het_has_no_note_key(self):
        data = TrioNode.get_zygosity_table_data()
        self.assertNotIn('note', data[TrioInheritance.COMPOUND_HET])

    # ── Cached label counts vs parent restriction ──────────────────────────────

    def test_only_compound_het_takes_a_parent(self):
        """Compound het is the only inheritance mode with an input; the rest are source nodes
        whose cohort-wide cached counts are valid (no parent restriction)."""
        self.assertTrue(self._make_node(TrioInheritance.COMPOUND_HET).has_input())
        for mode in (TrioInheritance.RECESSIVE, TrioInheritance.DENOVO, TrioInheritance.DOMINANT,
                     TrioInheritance.XLINKED_RECESSIVE, TrioInheritance.ALL_RECESSIVE,
                     TrioInheritance.ANY_AFFECTED):
            self.assertFalse(self._make_node(mode).has_input(), f"{mode} should be a source node")

    def test_compound_het_does_not_use_cohort_label_cache(self):
        """Compound het intersects its parent's queryset, so the cohort-wide stats cache (which
        ignores the parent) must NOT be consulted - otherwise it over-counts and trips the
        single-parent check in node_counts(). Regression for comp-het count > parent count."""
        node = self._make_node(TrioInheritance.COMPOUND_HET)
        with mock.patch(
                "analysis.models.nodes.sources.trio_node.get_cached_label_count_for_cohort",
                return_value=999999) as m:
            self.assertIsNone(node._get_cached_label_count(BuiltInFilters.TOTAL))
            m.assert_not_called()
