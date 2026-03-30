from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, TrioNode
from analysis.models.enums import TrioInheritance
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild, Variant
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
