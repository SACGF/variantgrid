from unittest import mock

from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, TrioNode
from analysis.models.enums import TrioInheritance
from analysis.tests.inheritance_node_mixin import InheritanceNodeTestsMixin, make_cohort_genotype
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import BuiltInFilters, GenomeBuild
from snpdb.models.models_cohort import CohortGenotypeCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestTrioNodeInheritance(InheritanceNodeTestsMixin, TestCase):
    """
    End-to-end tests for TrioNode inheritance filters.

    Creates real Variant + CohortGenotype records and verifies that the correct
    variants are returned (and wrong ones are excluded) for each inheritance mode.
    The modes shared with QuadNode come from InheritanceNodeTestsMixin; what's here is
    trio-specific - DOMINANT (the trio fixture is the one with an affected parent) and
    the compound-het count cache.

    Trio sample layout (from create_fake_trio → create_fake_cohort):
        packed_field_index 0 → proband
        packed_field_index 1 → mother  (mother_affected=True)
        packed_field_index 2 → father  (father_affected=False)

    samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
    """
    INHERITANCE = TrioInheritance
    NODE_CLASS = TrioNode
    REQUIRE_PARENT_ZYGOSITY = "require_zygosity"
    # samples_zygosity index: [proband, mother, father]
    SHARED_ZYGOSITIES = {
        "recessive_v": "OEE",                # proband=HOM_ALT, both parents=HET
        "denovo_v": "ERR",                   # proband=HET, both parents=HOM_REF
        "control_v": "RRR",                  # matches no inheritance mode
        "unknown_father_v": "OEU",           # recessive but father is a no-call
        "xlinked_v": "OER",                  # proband=HOM_ALT, mother=HET, on chrX
        "xlinked_unknown_mother_v": "OUR",   # XLR but mother is a no-call
    }

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

        cls.create_shared_variants(cls.cgc)

        # Dominant (mother affected): proband=HET, mother=HET, father=HOM_REF
        cls.dominant_v = cls._make_variant(3000, "EER")
        # Proband-only variant: proband=HET, mother=HOM_REF, father=HOM_REF
        cls.proband_only_v = cls._make_variant(6000, "ERR")
        # Mother-only variant: proband=HOM_REF, mother=HET, father=HOM_REF
        cls.mother_only_v = cls._make_variant(6100, "RER")
        # Father-only variant: proband=HOM_REF, mother=HOM_REF, father=HET
        cls.father_only_v = cls._make_variant(6200, "RRE")

    @classmethod
    def _make_variant(cls, position, samples_zygosity):
        variant = slowly_create_test_variant("3", position, "A", "T", cls.grch37)
        make_cohort_genotype(cls.cgc, variant, samples_zygosity)
        return variant

    def _make_node(self, inheritance, **kwargs):
        return TrioNode.objects.create(
            analysis=self.analysis, trio=self.trio,
            inheritance=inheritance, **kwargs
        )

    def _get_inheritance_errors(self, inheritance):
        return TrioNode.get_trio_inheritance_errors(self.trio, inheritance)

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

    def test_recessive_excludes_dominant_variant(self):
        node = self._make_node(TrioInheritance.RECESSIVE)
        self.assertNotIn(self.dominant_v.pk, self._filter_variants(node))

    def test_denovo_excludes_dominant_variant(self):
        node = self._make_node(TrioInheritance.DENOVO)
        self.assertNotIn(self.dominant_v.pk, self._filter_variants(node))

    def test_xlinked_excludes_dominant_variant(self):
        node = self._make_node(TrioInheritance.XLINKED_RECESSIVE)
        self.assertNotIn(self.dominant_v.pk, self._filter_variants(node))

    # ── Validation ────────────────────────────────────────────────────────────

    def test_dominant_with_affected_mother_no_errors(self):
        # trio has mother_affected=True
        self.assertEqual(self._get_inheritance_errors(TrioInheritance.DOMINANT), [])

    def test_dominant_no_affected_parent_raises_error(self):
        self.trio.mother_affected = False
        self.trio.father_affected = False
        errors = self._get_inheritance_errors(TrioInheritance.DOMINANT)
        self.assertGreater(len(errors), 0)
        self.assertIn("affected parent", errors[0].lower())
        self.trio.mother_affected = True  # restore
        self.trio.save()

    def test_xlinked_affected_mother_raises_error(self):
        # trio has mother_affected=True — xlinked requires unaffected mother
        errors = self._get_inheritance_errors(TrioInheritance.XLINKED_RECESSIVE)
        self.assertGreater(len(errors), 0)
        self.assertIn("mother", errors[0].lower())

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
                "analysis.models.nodes.sources.cohort_node.get_cached_label_count_for_cohort",
                return_value=999999) as m:
            self.assertIsNone(node._get_cached_label_count(BuiltInFilters.TOTAL))
            m.assert_not_called()
