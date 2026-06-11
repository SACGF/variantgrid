from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, QuadNode
from analysis.models.enums import QuadInheritance
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild, Variant
from snpdb.models.models_cohort import CohortGenotype, CohortGenotypeCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_quad
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestQuadNodeInheritance(TestCase):
    """
    End-to-end tests for QuadNode inheritance filters.

    Creates real Variant + CohortGenotype records and verifies that the correct
    variants are returned for each inheritance mode, including the sibling constraint
    that distinguishes Quad from Trio.

    Sample layout (from create_fake_quad):
        packed_field_index 0 → proband
        packed_field_index 1 → mother
        packed_field_index 2 → father
        packed_field_index 3 → sibling

    samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='testuser_quad')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        # ── Quad with unaffected sibling ─────────────────────────────────────
        cls.quad = create_fake_quad(user, cls.grch37, sibling_affected=False)
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cls.quad.cohort)

        # samples_zygosity index: [proband, mother, father, sibling]
        #
        # Recessive: proband=HOM_ALT, parents=HET, sibling=HOM_REF (carrier ok, not affected)
        cls.recessive_v = slowly_create_test_variant("3", 1000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc, cls.recessive_v, "OEER")

        # Denovo: proband=HET, parents=HOM_REF, sibling=HOM_REF
        cls.denovo_v = slowly_create_test_variant("3", 2000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc, cls.denovo_v, "ERRR")

        # X-linked: proband=HOM_ALT, mother=HET, father=any, sibling=HOM_REF — on chrX
        cls.xlinked_v = slowly_create_test_variant("X", 1000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc, cls.xlinked_v, "OERR")

        # Control: all HOM_REF — should match no mode
        cls.control_v = slowly_create_test_variant("3", 4000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc, cls.control_v, "RRRR")

        # Unknown-father autosomal variant: proband=HOM_ALT, mother=HET, father=UNKNOWN, sibling=HOM_REF
        cls.unknown_father_v = slowly_create_test_variant("3", 7000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc, cls.unknown_father_v, "OEUR")

        # X-linked with unknown mother: proband=HOM_ALT, mother=UNKNOWN, father=any, sibling=HOM_REF
        cls.xlinked_unknown_mother_v = slowly_create_test_variant("X", 2000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc, cls.xlinked_unknown_mother_v, "OURR")

        # ── Quad with affected sibling ────────────────────────────────────────
        cls.quad_aff = create_fake_quad(user, cls.grch37, sibling_affected=True)
        cls.cgc_aff = CohortGenotypeCollection.objects.get(cohort=cls.quad_aff.cohort)

        # Recessive (affected sibling): sibling must have HAS_VARIANT
        cls.recessive_aff_v = slowly_create_test_variant("3", 5000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc_aff, cls.recessive_aff_v, "OEEO")

        # Same pattern as unaffected-sibling recessive — sibling has HOM_REF, not HAS_VARIANT
        cls.recessive_sib_ref_v = slowly_create_test_variant("3", 6000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc_aff, cls.recessive_sib_ref_v, "OEER")

        # Proband-only variant (only proband has it): proband=HET, others HOM_REF
        cls.proband_only_v = slowly_create_test_variant("3", 8000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc_aff, cls.proband_only_v, "ERRR")

        # Sibling-only variant: proband=HOM_REF, mother=HOM_REF, father=HOM_REF, sibling=HET
        cls.sibling_only_v = slowly_create_test_variant("3", 8100, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc_aff, cls.sibling_only_v, "RRRE")

        # Mother-only variant: proband=HOM_REF, mother=HET, father=HOM_REF, sibling=HOM_REF
        cls.mother_only_v = slowly_create_test_variant("3", 8200, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc_aff, cls.mother_only_v, "RERR")

    @classmethod
    def _make_cg(cls, cgc, variant, samples_zygosity):
        n = len(samples_zygosity)
        CohortGenotype.objects.create(
            collection=cgc,
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

    def _make_node(self, inheritance, quad=None, **kwargs):
        return QuadNode.objects.create(
            analysis=self.analysis,
            quad=quad or self.quad,
            inheritance=inheritance,
            **kwargs
        )

    def _filter_variants(self, node):
        """Apply node's inheritance filter against real Variant data; return matching PKs."""
        arg_q_dict = node._get_node_arg_q_dict()
        cgc = node.quad.cohort.cohort_genotype_collection
        qs = Variant.objects.annotate(**cgc.get_annotation_kwargs())
        alias = cgc.cohortgenotype_alias
        if alias in arg_q_dict:
            for q in arg_q_dict[alias].values():
                qs = qs.filter(q)
        if None in arg_q_dict:
            for q in arg_q_dict[None].values():
                qs = qs.filter(q)
        return set(qs.values_list('pk', flat=True))

    # ── Recessive (unaffected sibling) ────────────────────────────────────────

    def test_recessive_matches_recessive_variant(self):
        node = self._make_node(QuadInheritance.RECESSIVE)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_recessive_excludes_other_variants(self):
        node = self._make_node(QuadInheritance.RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── Denovo ────────────────────────────────────────────────────────────────

    def test_denovo_matches_denovo_variant(self):
        node = self._make_node(QuadInheritance.DENOVO)
        self.assertIn(self.denovo_v.pk, self._filter_variants(node))

    def test_denovo_excludes_other_variants(self):
        node = self._make_node(QuadInheritance.DENOVO)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── X-linked recessive ────────────────────────────────────────────────────

    def test_xlinked_matches_x_chromosome_variant(self):
        node = self._make_node(QuadInheritance.XLINKED_RECESSIVE)
        self.assertIn(self.xlinked_v.pk, self._filter_variants(node))

    def test_xlinked_excludes_autosomal_variants(self):
        node = self._make_node(QuadInheritance.XLINKED_RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.denovo_v.pk, ids)

    # ── Sibling constraint (key distinction vs Trio) ──────────────────────────

    def test_affected_sibling_recessive_matches_sibling_with_variant(self):
        """Affected-sibling recessive requires sibling to also carry the variant."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad_aff)
        ids = self._filter_variants(node)
        self.assertIn(self.recessive_aff_v.pk, ids)

    def test_affected_sibling_recessive_excludes_sibling_without_variant(self):
        """With sibling_affected=True, a variant where sibling=HOM_REF is excluded."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad_aff)
        ids = self._filter_variants(node)
        self.assertNotIn(self.recessive_sib_ref_v.pk, ids)

    def test_unaffected_sibling_recessive_includes_sibling_as_carrier(self):
        """With sibling_affected=False, sibling can be a carrier (HOM_REF or HET)."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad)
        ids = self._filter_variants(node)
        self.assertIn(self.recessive_v.pk, ids)  # sibling='R' (HOM_REF) is allowed

    # ── Validation ────────────────────────────────────────────────────────────

    def test_dominant_no_affected_parent_raises_error(self):
        # Both quads have mother_affected=False, father_affected=False
        errors = QuadNode.get_quad_inheritance_errors(self.quad, QuadInheritance.DOMINANT)
        self.assertGreater(len(errors), 0)

    def test_xlinked_affected_mother_raises_error(self):
        self.quad.mother_affected = True
        self.quad.save()
        errors = QuadNode.get_quad_inheritance_errors(self.quad, QuadInheritance.XLINKED_RECESSIVE)
        self.assertGreater(len(errors), 0)
        self.quad.mother_affected = False
        self.quad.save()

    # ── CompHet node type ─────────────────────────────────────────────────────

    def test_compound_het_requires_one_parent_input(self):
        node = self._make_node(QuadInheritance.COMPOUND_HET)
        self.assertEqual(node.min_inputs, 1)
        self.assertEqual(node.max_inputs, 1)

    def test_quality_filter_ignores_foreign_ancestor_samples(self):
        """Regression: COMPOUND_HET is the only mode with a parent input, so get_sample_ids()
        includes ancestor samples that may not belong to this node's cohort. Per-sample quality
        filters must skip those rather than raising KeyError in get_array_index_for_sample_id."""
        node = self._make_node(QuadInheritance.COMPOUND_HET, min_dp=30)

        quad_sample_ids = list(self.quad.cohort.get_sample_ids())
        # A sample from the *other* quad's cohort - simulates a parent node feeding a foreign sample
        foreign_sample_id = self.quad_aff.proband.sample_id
        self.assertNotIn(foreign_sample_id, quad_sample_ids)

        with patch.object(node, "get_sample_ids", return_value=quad_sample_ids + [foreign_sample_id]):
            cohort, arg_q_dict = node.get_cohort_and_arg_q_dict()  # Must not raise KeyError

        self.assertEqual(cohort, self.quad.cohort)
        # min_dp filter is still applied (to the quad's own samples)
        self.assertIn(self.cgc.cohortgenotype_alias, arg_q_dict)

    def test_simple_modes_are_source_nodes(self):
        for mode in [QuadInheritance.RECESSIVE, QuadInheritance.DENOVO,
                     QuadInheritance.DOMINANT, QuadInheritance.XLINKED_RECESSIVE]:
            node = self._make_node(mode)
            self.assertEqual(node.max_inputs, 0, f"{mode} should be a source node")

    # ── Clone ─────────────────────────────────────────────────────────────────

    def test_clone_matches_same_variants(self):
        node = self._make_node(QuadInheritance.RECESSIVE)
        clone = node.save_clone()
        self.assertEqual(self._filter_variants(node), self._filter_variants(clone))

    # ── ALL_RECESSIVE (AR ∪ XLR) ──────────────────────────────────────────────

    def test_all_recessive_matches_autosomal_recessive_variant(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_all_recessive_matches_xlinked_variant(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE)
        self.assertIn(self.xlinked_v.pk, self._filter_variants(node))

    def test_all_recessive_excludes_unrelated_variants(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE)
        ids = self._filter_variants(node)
        self.assertNotIn(self.denovo_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    def test_all_recessive_require_zygosity_excludes_unknown_father_on_autosome(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE, require_zygosity=True)
        self.assertNotIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_all_recessive_no_require_zygosity_includes_unknown_father_on_autosome(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE, require_zygosity=False)
        self.assertIn(self.unknown_father_v.pk, self._filter_variants(node))

    def test_all_recessive_require_zygosity_excludes_unknown_mother_on_xlr(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE, require_zygosity=True)
        self.assertNotIn(self.xlinked_unknown_mother_v.pk, self._filter_variants(node))

    def test_all_recessive_no_require_zygosity_includes_unknown_mother_on_xlr(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE, require_zygosity=False)
        self.assertIn(self.xlinked_unknown_mother_v.pk, self._filter_variants(node))

    def test_all_recessive_is_source_node(self):
        node = self._make_node(QuadInheritance.ALL_RECESSIVE)
        self.assertEqual(node.max_inputs, 0)

    def test_zygosity_table_all_recessive_two_line_cells(self):
        data = QuadNode.get_zygosity_table_data()
        self.assertIn('AR:', data[QuadInheritance.ALL_RECESSIVE]['mother'])
        self.assertIn('XLR:', data[QuadInheritance.ALL_RECESSIVE]['mother'])
        self.assertIn('AR:', data[QuadInheritance.ALL_RECESSIVE]['sibling'])

    def test_zygosity_table_all_recessive_other_filters_mentions_chr_x(self):
        data = QuadNode.get_zygosity_table_data()
        self.assertIn('Chr X', data[QuadInheritance.ALL_RECESSIVE]['other_filters_mother'])

    # ── ANY_AFFECTED ──────────────────────────────────────────────────────────

    def test_any_affected_sibling_affected_includes_proband_only_variant(self):
        # quad_aff has sibling_affected=True
        node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad_aff)
        self.assertIn(self.proband_only_v.pk, self._filter_variants(node))

    def test_any_affected_sibling_affected_includes_sibling_only_variant(self):
        node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad_aff)
        self.assertIn(self.sibling_only_v.pk, self._filter_variants(node))

    def test_any_affected_sibling_affected_excludes_mother_only_variant(self):
        # mother is unaffected on quad_aff
        node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad_aff)
        self.assertNotIn(self.mother_only_v.pk, self._filter_variants(node))

    def test_any_affected_sibling_unaffected_excludes_sibling_only_variant(self):
        # cls.quad has sibling_affected=False
        # sibling_only_v lives in cgc_aff; create one in the unaffected quad's cgc
        sib_only_unaff_v = slowly_create_test_variant("3", 9000, "A", "T", self.grch37)
        self._make_cg(self.cgc, sib_only_unaff_v, "RRRE")
        node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad)
        self.assertNotIn(sib_only_unaff_v.pk, self._filter_variants(node))

    def test_any_affected_sibling_unaffected_includes_proband_only_variant(self):
        proband_only_unaff_v = slowly_create_test_variant("3", 9100, "A", "T", self.grch37)
        self._make_cg(self.cgc, proband_only_unaff_v, "ERRR")
        node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad)
        self.assertIn(proband_only_unaff_v.pk, self._filter_variants(node))

    def test_any_affected_with_affected_mother_includes_mother_only(self):
        self.quad_aff.mother_affected = True
        self.quad_aff.save()
        try:
            node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad_aff)
            ids = self._filter_variants(node)
            self.assertIn(self.mother_only_v.pk, ids)
        finally:
            self.quad_aff.mother_affected = False
            self.quad_aff.save()

    def test_any_affected_is_source_node(self):
        node = self._make_node(QuadInheritance.ANY_AFFECTED)
        self.assertEqual(node.max_inputs, 0)

    def test_any_affected_always_valid_no_errors(self):
        errors = QuadNode.get_quad_inheritance_errors(self.quad, QuadInheritance.ANY_AFFECTED)
        self.assertEqual(errors, [])

    def test_zygosity_table_any_affected_proband_has_variant(self):
        data = QuadNode.get_zygosity_table_data()
        proband_value = data[QuadInheritance.ANY_AFFECTED]['proband']
        self.assertTrue(proband_value)

    # ── Zygosity table "Other Filters" column ─────────────────────────────────

    def test_zygosity_table_xlinked_has_chr_x_other_filter(self):
        data = QuadNode.get_zygosity_table_data()
        self.assertEqual(data[QuadInheritance.XLINKED_RECESSIVE]['other_filters_mother'], "Chr X only")

    def test_zygosity_table_recessive_has_no_other_filter_keys(self):
        data = QuadNode.get_zygosity_table_data()
        keys = [k for k in data[QuadInheritance.RECESSIVE] if k.startswith('other_filters_')]
        self.assertEqual(keys, [])

    def test_zygosity_table_compound_het_has_gene_constraint_in_other_filters(self):
        data = QuadNode.get_zygosity_table_data()
        self.assertIn("gene", data[QuadInheritance.COMPOUND_HET]['other_filters_mother'])

    def test_zygosity_table_compound_het_has_no_note_key(self):
        data = QuadNode.get_zygosity_table_data()
        self.assertNotIn('note', data[QuadInheritance.COMPOUND_HET])
