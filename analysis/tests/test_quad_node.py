from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, QuadNode
from analysis.models.enums import QuadInheritance
from analysis.tests.inheritance_node_mixin import InheritanceNodeTestsMixin, make_cohort_genotype
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild
from snpdb.models.models_cohort import CohortGenotypeCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_quad
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestQuadNodeInheritance(InheritanceNodeTestsMixin, TestCase):
    """
    End-to-end tests for QuadNode inheritance filters.

    The inheritance modes shared with TrioNode come from InheritanceNodeTestsMixin (run here
    over the 4-sample layout); what's here is the sibling constraint that distinguishes Quad
    from Trio, plus the compound-het foreign-sample regression.

    Sample layout (from create_fake_quad) - both parents are unaffected:
        packed_field_index 0 → proband
        packed_field_index 1 → mother
        packed_field_index 2 → father
        packed_field_index 3 → sibling

    samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
    """
    INHERITANCE = QuadInheritance
    NODE_CLASS = QuadNode
    REQUIRE_PARENT_ZYGOSITY = "require_parent_zygosity"
    # samples_zygosity index: [proband, mother, father, sibling]
    SHARED_ZYGOSITIES = {
        "recessive_v": "OEER",               # sibling=HOM_REF - a carrier is ok, not affected
        "denovo_v": "ERRR",
        "control_v": "RRRR",
        "unknown_father_v": "OEUR",
        "xlinked_v": "OERR",
        "xlinked_unknown_mother_v": "OURR",
    }

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
        cls.create_shared_variants(cls.cgc)

        # ── Quad with affected sibling ────────────────────────────────────────
        cls.quad_aff = create_fake_quad(user, cls.grch37, sibling_affected=True)
        cls.cgc_aff = CohortGenotypeCollection.objects.get(cohort=cls.quad_aff.cohort)

        # Recessive (affected sibling): sibling must have HAS_VARIANT
        cls.recessive_aff_v = cls._make_variant(cls.cgc_aff, 5000, "OEEO")
        # Same pattern as unaffected-sibling recessive — sibling has HOM_REF, not HAS_VARIANT
        cls.recessive_sib_ref_v = cls._make_variant(cls.cgc_aff, 6000, "OEER")
        # Proband-only variant: proband=HET, others HOM_REF
        cls.proband_only_v = cls._make_variant(cls.cgc_aff, 8000, "ERRR")
        # Sibling-only variant: sibling=HET, everyone else HOM_REF
        cls.sibling_only_v = cls._make_variant(cls.cgc_aff, 8100, "RRRE")
        # Mother-only variant: mother=HET, everyone else HOM_REF
        cls.mother_only_v = cls._make_variant(cls.cgc_aff, 8200, "RERR")
        # Recessive w/ unknown sibling: proband=HOM_ALT, parents=HET, sibling=UNKNOWN.
        # On the affected-sibling quad the sibling must HAS_VARIANT, so this only passes
        # when sibling zygosity is not required.
        cls.recessive_sib_unknown_v = cls._make_variant(cls.cgc_aff, 8300, "OEEU")

    @classmethod
    def _make_variant(cls, cgc, position, samples_zygosity):
        variant = slowly_create_test_variant("3", position, "A", "T", cls.grch37)
        make_cohort_genotype(cgc, variant, samples_zygosity)
        return variant

    def _make_node(self, inheritance, quad=None, **kwargs):
        return QuadNode.objects.create(
            analysis=self.analysis,
            quad=quad or self.quad,
            inheritance=inheritance,
            **kwargs
        )

    def _get_inheritance_errors(self, inheritance):
        return QuadNode.get_quad_inheritance_errors(self.quad, inheritance)

    # ── Sibling constraint (key distinction vs Trio) ──────────────────────────

    def test_affected_sibling_recessive_matches_sibling_with_variant(self):
        """Affected-sibling recessive requires sibling to also carry the variant."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad_aff)
        self.assertIn(self.recessive_aff_v.pk, self._filter_variants(node))

    def test_affected_sibling_recessive_excludes_sibling_without_variant(self):
        """With sibling_affected=True, a variant where sibling=HOM_REF is excluded."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad_aff)
        self.assertNotIn(self.recessive_sib_ref_v.pk, self._filter_variants(node))

    def test_unaffected_sibling_recessive_includes_sibling_as_carrier(self):
        """With sibling_affected=False, sibling can be a carrier (HOM_REF or HET)."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))  # sibling='R' is allowed

    # ── require_sibling_zygosity ──────────────────────────────────────────────

    def test_require_sibling_zygosity_true_excludes_unknown_sibling(self):
        """Default (require_sibling_zygosity=True) excludes a sibling with no genotype call."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad_aff,
                               require_sibling_zygosity=True)
        self.assertNotIn(self.recessive_sib_unknown_v.pk, self._filter_variants(node))

    def test_require_sibling_zygosity_false_includes_unknown_sibling(self):
        """With require_sibling_zygosity=False, a no-call sibling is allowed through."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad_aff,
                               require_sibling_zygosity=False)
        self.assertIn(self.recessive_sib_unknown_v.pk, self._filter_variants(node))

    def test_require_sibling_zygosity_independent_of_parent_toggle(self):
        """Requiring sibling zygosity but not parent zygosity still excludes a no-call sibling."""
        node = self._make_node(QuadInheritance.RECESSIVE, quad=self.quad_aff,
                               require_parent_zygosity=False, require_sibling_zygosity=True)
        self.assertNotIn(self.recessive_sib_unknown_v.pk, self._filter_variants(node))

    # ── Validation ────────────────────────────────────────────────────────────

    def test_dominant_no_affected_parent_raises_error(self):
        # Both quads have mother_affected=False, father_affected=False
        errors = self._get_inheritance_errors(QuadInheritance.DOMINANT)
        self.assertGreater(len(errors), 0)

    def test_xlinked_affected_mother_raises_error(self):
        self.quad.mother_affected = True
        self.quad.save()
        try:
            errors = self._get_inheritance_errors(QuadInheritance.XLINKED_RECESSIVE)
            self.assertGreater(len(errors), 0)
        finally:
            self.quad.mother_affected = False
            self.quad.save()

    # ── CompHet node type ─────────────────────────────────────────────────────

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
        sib_only_unaff_v = self._make_variant(self.cgc, 9000, "RRRE")
        node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad)
        self.assertNotIn(sib_only_unaff_v.pk, self._filter_variants(node))

    def test_any_affected_sibling_unaffected_includes_proband_only_variant(self):
        proband_only_unaff_v = self._make_variant(self.cgc, 9100, "ERRR")
        node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad)
        self.assertIn(proband_only_unaff_v.pk, self._filter_variants(node))

    def test_any_affected_with_affected_mother_includes_mother_only(self):
        self.quad_aff.mother_affected = True
        self.quad_aff.save()
        try:
            node = self._make_node(QuadInheritance.ANY_AFFECTED, quad=self.quad_aff)
            self.assertIn(self.mother_only_v.pk, self._filter_variants(node))
        finally:
            self.quad_aff.mother_affected = False
            self.quad_aff.save()

    # ── Zygosity table ────────────────────────────────────────────────────────

    def test_zygosity_table_all_recessive_sibling_two_line_cell(self):
        data = QuadNode.get_zygosity_table_data()
        self.assertIn('AR:', data[QuadInheritance.ALL_RECESSIVE]['sibling'])
