"""
End-to-end tests for DuoNode inheritance filters.

DuoNode doesn't share InheritanceNodeTestsMixin with Trio/Quad: with one parent the modes are a
different set (there's no Denovo, only "absent in parent") and the zygosity table is keyed on
'parent' rather than mother/father.

Sample layout (from create_fake_duo):
    packed_field_index 0 → proband
    packed_field_index 1 → parent

samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, DuoNode
from analysis.models.enums import DuoInheritance
from analysis.tests.inheritance_node_mixin import make_cohort_genotype
from annotation.fake_annotation import get_fake_annotation_version
from patients.models_enums import Sex
from snpdb.models import Duo, DuoRelationship, GenomeBuild, Variant
from snpdb.models.models_cohort import CohortGenotypeCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_duo
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestDuoNodeInheritance(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='testuser_duo')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        cls.duo = create_fake_duo(user, cls.grch37, relationship=DuoRelationship.MOTHER)
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cls.duo.cohort)
        # The same two samples read as a father duo - the modes that care take the relationship
        cls.duo_father = Duo.objects.create(name="test_duo_father", user=user, cohort=cls.duo.cohort,
                                            proband=cls.duo.proband, parent=cls.duo.parent,
                                            relationship=DuoRelationship.FATHER)

        # samples_zygosity index: [proband, parent]
        cls.recessive_v = cls._make_variant("3", 1000, "OE")       # proband HOM_ALT, parent HET
        cls.absent_v = cls._make_variant("3", 2000, "ER")          # proband HET, parent HOM_REF
        cls.control_v = cls._make_variant("3", 4000, "RR")
        cls.unknown_parent_v = cls._make_variant("3", 5000, "OU")  # recessive but parent is a no-call
        cls.parent_only_v = cls._make_variant("3", 6000, "RE")
        cls.dominant_v = cls._make_variant("3", 7000, "EE")        # both have it
        cls.xlinked_v = cls._make_variant("X", 1000, "OE")

    @classmethod
    def _make_variant(cls, chrom, position, samples_zygosity):
        variant = slowly_create_test_variant(chrom, position, "A", "T", cls.grch37)
        make_cohort_genotype(cls.cgc, variant, samples_zygosity)
        return variant

    def _make_node(self, inheritance, duo=None, **kwargs):
        return DuoNode.objects.create(analysis=self.analysis, duo=duo or self.duo,
                                      inheritance=inheritance, **kwargs)

    def _filter_variants(self, node):
        """Apply node's inheritance filter against real Variant data; return matching PKs."""
        arg_q_dict = node._get_node_arg_q_dict()
        cgc = node._get_cohort().cohort_genotype_collection
        qs = Variant.objects.annotate(**cgc.get_annotation_kwargs())
        for alias in (cgc.cohortgenotype_alias, None):
            for q in arg_q_dict.get(alias, {}).values():
                qs = qs.filter(q)
        return set(qs.values_list('pk', flat=True))

    # ── Recessive ─────────────────────────────────────────────────────────────

    def test_recessive_matches_recessive_variant(self):
        node = self._make_node(DuoInheritance.RECESSIVE)
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_recessive_excludes_other_variants(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.RECESSIVE))
        self.assertNotIn(self.absent_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── Absent in parent ──────────────────────────────────────────────────────

    def test_absent_in_parent_matches_variant_the_parent_lacks(self):
        node = self._make_node(DuoInheritance.ABSENT_IN_PARENT)
        self.assertIn(self.absent_v.pk, self._filter_variants(node))

    def test_absent_in_parent_excludes_inherited_variants(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.ABSENT_IN_PARENT))
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.dominant_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    def test_absent_in_parent_always_warns_about_the_missing_parent(self):
        node = self._make_node(DuoInheritance.ABSENT_IN_PARENT)
        warnings = node.get_warnings()
        self.assertTrue(warnings)
        self.assertIn("father", warnings[0].lower())  # the duo has the mother

    def test_absent_in_parent_warning_names_the_missing_mother(self):
        node = self._make_node(DuoInheritance.ABSENT_IN_PARENT, duo=self.duo_father)
        self.assertIn("mother", node.get_warnings()[0].lower())

    def test_other_modes_do_not_warn(self):
        self.assertEqual(self._make_node(DuoInheritance.RECESSIVE).get_warnings(), [])

    # ── Dominant ──────────────────────────────────────────────────────────────

    def test_dominant_unaffected_parent_requires_variant_absent_in_parent(self):
        node = self._make_node(DuoInheritance.DOMINANT)
        ids = self._filter_variants(node)
        self.assertIn(self.absent_v.pk, ids)
        self.assertNotIn(self.dominant_v.pk, ids)

    def test_dominant_affected_parent_matches_shared_variant(self):
        self.duo.parent_affected = True
        try:
            node = self._make_node(DuoInheritance.DOMINANT)
            ids = self._filter_variants(node)
            self.assertIn(self.dominant_v.pk, ids)
            self.assertNotIn(self.absent_v.pk, ids)
        finally:
            self.duo.parent_affected = False

    # ── X-linked recessive ────────────────────────────────────────────────────

    def test_xlinked_matches_x_chromosome_variant(self):
        node = self._make_node(DuoInheritance.XLINKED_RECESSIVE)
        self.assertIn(self.xlinked_v.pk, self._filter_variants(node))

    def test_xlinked_excludes_autosomal_variants(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.XLINKED_RECESSIVE))
        self.assertNotIn(self.recessive_v.pk, ids)
        self.assertNotIn(self.absent_v.pk, ids)

    # ── All recessive (AR ∪ XLR) ──────────────────────────────────────────────

    def test_all_recessive_with_mother_matches_both_branches(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.ALL_RECESSIVE))
        self.assertIn(self.recessive_v.pk, ids)
        self.assertIn(self.xlinked_v.pk, ids)

    def test_all_recessive_with_father_drops_the_xlinked_branch(self):
        """ chrX doesn't come from the father, so there's no XLR arm to OR in """
        node = self._make_node(DuoInheritance.ALL_RECESSIVE, duo=self.duo_father)
        self.assertNotIn("XLR", node._get_method_summary())
        self.assertIn(self.recessive_v.pk, self._filter_variants(node))

    def test_all_recessive_with_mother_keeps_the_xlinked_branch(self):
        node = self._make_node(DuoInheritance.ALL_RECESSIVE)
        self.assertIn("XLR", node._get_method_summary())

    def test_all_recessive_excludes_unrelated_variants(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.ALL_RECESSIVE))
        self.assertNotIn(self.absent_v.pk, ids)
        self.assertNotIn(self.control_v.pk, ids)

    # ── Any affected ──────────────────────────────────────────────────────────

    def test_any_affected_unaffected_parent_collapses_to_proband(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.ANY_AFFECTED))
        self.assertIn(self.absent_v.pk, ids)
        self.assertNotIn(self.parent_only_v.pk, ids)

    def test_any_affected_affected_parent_includes_parent_only_variant(self):
        self.duo.parent_affected = True
        try:
            ids = self._filter_variants(self._make_node(DuoInheritance.ANY_AFFECTED))
            self.assertIn(self.parent_only_v.pk, ids)
        finally:
            self.duo.parent_affected = False

    # ── require_zygosity (parent only) ────────────────────────────────────────

    def test_require_zygosity_true_excludes_parent_no_call(self):
        node = self._make_node(DuoInheritance.RECESSIVE, require_zygosity=True)
        self.assertNotIn(self.unknown_parent_v.pk, self._filter_variants(node))

    def test_require_zygosity_false_allows_parent_no_call(self):
        node = self._make_node(DuoInheritance.RECESSIVE, require_zygosity=False)
        self.assertIn(self.unknown_parent_v.pk, self._filter_variants(node))

    def test_require_zygosity_false_still_requires_proband_zygosity(self):
        """ #947 - widening applies to the parent only, so a no-call proband is still excluded """
        proband_no_call_v = self._make_variant("3", 9000, "UE")
        node = self._make_node(DuoInheritance.RECESSIVE, require_zygosity=False)
        self.assertNotIn(proband_no_call_v.pk, self._filter_variants(node))

    # ── Node type ─────────────────────────────────────────────────────────────

    def test_compound_het_requires_one_parent_input(self):
        node = self._make_node(DuoInheritance.COMPOUND_HET)
        self.assertEqual(node.min_inputs, 1)
        self.assertEqual(node.max_inputs, 1)

    def test_other_modes_are_source_nodes(self):
        for mode in DuoInheritance:
            if mode != DuoInheritance.COMPOUND_HET:
                self.assertEqual(self._make_node(mode).max_inputs, 0, f"{mode} should be a source node")

    # ── Validation ────────────────────────────────────────────────────────────

    def test_dominant_unaffected_parent_raises_error(self):
        errors = DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.DOMINANT)
        self.assertGreater(len(errors), 0)
        self.assertIn("affected parent", errors[0].lower())

    def test_dominant_affected_parent_no_errors(self):
        self.duo.parent_affected = True
        try:
            self.assertEqual(DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.DOMINANT), [])
        finally:
            self.duo.parent_affected = False

    def test_xlinked_with_father_raises_error(self):
        errors = DuoNode.get_duo_inheritance_errors(self.duo_father, DuoInheritance.XLINKED_RECESSIVE)
        self.assertGreater(len(errors), 0)
        self.assertIn("mother", errors[0].lower())

    def test_xlinked_with_mother_no_errors(self):
        self.assertEqual(DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.XLINKED_RECESSIVE), [])

    def test_xlinked_affected_mother_raises_error(self):
        self.duo.parent_affected = True
        try:
            errors = DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.XLINKED_RECESSIVE)
            self.assertGreater(len(errors), 0)
            self.assertIn("unaffected mother", errors[0].lower())
        finally:
            self.duo.parent_affected = False

    def test_xlinked_female_proband_raises_error(self):
        self.duo.proband_sex = Sex.FEMALE
        try:
            errors = DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.XLINKED_RECESSIVE)
            self.assertGreater(len(errors), 0)
            self.assertIn("female proband", errors[0].lower())
        finally:
            self.duo.proband_sex = None

    def test_any_affected_always_valid_no_errors(self):
        self.assertEqual(DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.ANY_AFFECTED), [])

    # ── Zygosity table ────────────────────────────────────────────────────────

    def test_zygosity_table_all_recessive_mother_has_xlinked_branch(self):
        data = DuoNode.get_zygosity_table_data()
        self.assertIn('XLR:', data[DuoInheritance.ALL_RECESSIVE][f'parent_{DuoRelationship.MOTHER}'])

    def test_zygosity_table_all_recessive_father_has_no_xlinked_branch(self):
        data = DuoNode.get_zygosity_table_data()
        self.assertNotIn('XLR:', data[DuoInheritance.ALL_RECESSIVE][f'parent_{DuoRelationship.FATHER}'])

    def test_zygosity_table_xlinked_has_chr_x_other_filter(self):
        data = DuoNode.get_zygosity_table_data()
        self.assertEqual(data[DuoInheritance.XLINKED_RECESSIVE]['other_filters_parent'], "Chr X only")

    def test_zygosity_table_compound_het_says_one_hit_from_the_parent(self):
        data = DuoNode.get_zygosity_table_data()
        self.assertIn("gene", data[DuoInheritance.COMPOUND_HET]['other_filters_parent'])

    def test_zygosity_table_any_affected_parent_only_counted_when_affected(self):
        data = DuoNode.get_zygosity_table_data()
        entry = data[DuoInheritance.ANY_AFFECTED]
        self.assertTrue(entry['parent_affected'])
        self.assertNotEqual(entry['parent_affected'], entry['parent_unaffected'])
