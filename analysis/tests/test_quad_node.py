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

        # ── Quad with affected sibling ────────────────────────────────────────
        cls.quad_aff = create_fake_quad(user, cls.grch37, sibling_affected=True)
        cls.cgc_aff = CohortGenotypeCollection.objects.get(cohort=cls.quad_aff.cohort)

        # Recessive (affected sibling): sibling must have HAS_VARIANT
        cls.recessive_aff_v = slowly_create_test_variant("3", 5000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc_aff, cls.recessive_aff_v, "OEEO")

        # Same pattern as unaffected-sibling recessive — sibling has HOM_REF, not HAS_VARIANT
        cls.recessive_sib_ref_v = slowly_create_test_variant("3", 6000, "A", "T", cls.grch37)
        cls._make_cg(cls.cgc_aff, cls.recessive_sib_ref_v, "OEER")

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
