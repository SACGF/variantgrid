"""
End-to-end tests for DuoNode inheritance filters.

DuoNode doesn't share InheritanceNodeTestsMixin with Trio/Quad: with one relative the modes are a
different set (there's no Denovo, only "absent in parent") and the zygosity table is keyed on
'relative' rather than mother/father.

Sample layout (from create_fake_duo):
    packed_field_index 0 → proband
    packed_field_index 1 → relative

samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT, U=UNKNOWN, .=MISSING
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis, DuoNode
from analysis.models.enums import DuoInheritance
from analysis.models.nodes.family_inheritance import MOSAIC_JOINT_CALL_WARNING
from analysis.models.nodes.sources.duo_node import COMP_HET_SIBLING_UNPHASED
from analysis.tests.inheritance_node_mixin import DEFAULT_GENOTYPE_VALUES, make_cohort_genotype
from annotation.fake_annotation import get_fake_annotation_version
from patients.models_enums import Sex, Zygosity
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
        # The same two samples read as other kinds of duo - the modes that care take the relationship
        cls.duo_father = cls._make_duo(user, "test_duo_father", DuoRelationship.FATHER)
        cls.duo_sibling = cls._make_duo(user, "test_duo_sibling", DuoRelationship.SIBLING)
        cls.duo_sibling_affected = cls._make_duo(user, "test_duo_sibling_affected",
                                                 DuoRelationship.SIBLING, relative_affected=True)

        # samples_zygosity index: [proband, relative]
        cls.recessive_v = cls._make_variant("3", 1000, "OE")       # proband HOM_ALT, relative HET
        cls.absent_v = cls._make_variant("3", 2000, "ER")          # proband HET, relative HOM_REF
        cls.control_v = cls._make_variant("3", 4000, "RR")
        cls.unknown_parent_v = cls._make_variant("3", 5000, "OU")  # recessive but relative is a no-call
        cls.parent_only_v = cls._make_variant("3", 6000, "RE")
        cls.dominant_v = cls._make_variant("3", 7000, "EE")        # both have it
        cls.xlinked_v = cls._make_variant("X", 1000, "OE")
        cls.both_hom_alt_v = cls._make_variant("3", 10000, "OO")   # both homozygous
        cls.x_both_hom_alt_v = cls._make_variant("X", 10000, "OO")

        # Mosaic parent (#1830) - AD/AF matter as much as the call. [proband, parent]
        # Parent called HOM_REF but with 3 alt reads at 6%
        cls.mosaic_parent_v = cls._make_variant("3", 8000, "ER", [15, 3], [0.5, 0.06])
        # Parent called HET, but at 15% - too low to be constitutional
        cls.mosaic_het_parent_v = cls._make_variant("3", 8100, "EE", [15, 4], [0.5, 0.15])
        # Parent is a full HET - inherited dominant, not mosaic
        cls.inherited_het_v = cls._make_variant("3", 8200, "EE", [15, 14], [0.5, 0.48])
        # A single alt read in the parent - below the evidence threshold
        cls.single_alt_read_v = cls._make_variant("3", 8300, "ER", [15, 1], [0.5, 0.02])
        # Mosaic parent, but the proband doesn't carry it
        cls.parent_mosaic_no_proband_v = cls._make_variant("3", 8400, "RR", [0, 3], [0.0, 0.06])
        # AF FORMAT field present but no value for this record - the missing number
        cls.mosaic_missing_af_v = cls._make_variant("3", 8500, "ER", [15, 3], [-1, -1])
        # VCF had no AF FORMAT field at all, so the whole column is NULL
        cls.mosaic_null_af_v = cls._make_variant("3", 8600, "ER", [15, 3], allele_frequency=None)

    @classmethod
    def _make_duo(cls, user, name, relationship, relative_affected=False) -> Duo:
        return Duo.objects.create(name=name, user=user, cohort=cls.duo.cohort, proband=cls.duo.proband,
                                  relative=cls.duo.relative, relationship=relationship,
                                  relative_affected=relative_affected)

    @classmethod
    def _make_variant(cls, chrom, position, samples_zygosity, allele_depth=DEFAULT_GENOTYPE_VALUES,
                      allele_frequency=DEFAULT_GENOTYPE_VALUES):
        variant = slowly_create_test_variant(chrom, position, "A", "T", cls.grch37)
        make_cohort_genotype(cls.cgc, variant, samples_zygosity, allele_depth, allele_frequency)
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
        self.duo.relative_affected = True
        try:
            node = self._make_node(DuoInheritance.DOMINANT)
            ids = self._filter_variants(node)
            self.assertIn(self.dominant_v.pk, ids)
            self.assertNotIn(self.absent_v.pk, ids)
        finally:
            self.duo.relative_affected = False

    # ── Dominant (mosaic parent) ──────────────────────────────────────────────

    def test_mosaic_matches_alt_reads_in_a_hom_ref_called_parent(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        self.assertIn(self.mosaic_parent_v.pk, self._filter_variants(node))

    def test_mosaic_matches_low_vaf_het_called_parent(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        self.assertIn(self.mosaic_het_parent_v.pk, self._filter_variants(node))

    def test_mosaic_excludes_full_het_parent(self):
        """ A parent at 48% is a constitutional het - that's plain dominant, not mosaic """
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        self.assertNotIn(self.inherited_het_v.pk, self._filter_variants(node))

    def test_mosaic_excludes_parent_below_min_alt_reads(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        self.assertNotIn(self.single_alt_read_v.pk, self._filter_variants(node))

    def test_mosaic_excludes_variants_the_proband_lacks(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        self.assertNotIn(self.parent_mosaic_no_proband_v.pk, self._filter_variants(node))

    def test_mosaic_matches_when_the_allele_frequency_is_missing(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        self.assertIn(self.mosaic_missing_af_v.pk, self._filter_variants(node))

    def test_mosaic_matches_when_the_vcf_has_no_allele_frequency_field(self):
        """ AD carries the mode on its own - AF is only stored when the VCF has the FORMAT field """
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        self.assertIn(self.mosaic_null_af_v.pk, self._filter_variants(node))

    def test_mosaic_min_alt_reads_raises_the_bar(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT, mosaic_min_alt_reads=4)
        self.assertNotIn(self.mosaic_parent_v.pk, self._filter_variants(node))

    def test_mosaic_max_af_lets_a_higher_vaf_parent_in(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT, mosaic_max_af=0.6)
        self.assertIn(self.inherited_het_v.pk, self._filter_variants(node))

    def test_mosaic_joint_called_cohort_skips_the_joint_call_warning(self):
        warnings = self._make_node(DuoInheritance.MOSAIC_PARENT).get_warnings()
        self.assertTrue(warnings)
        self.assertNotIn(MOSAIC_JOINT_CALL_WARNING, warnings)

    def test_mosaic_multi_vcf_cohort_warns_it_needs_a_joint_call(self):
        node = self._make_node(DuoInheritance.MOSAIC_PARENT)
        node.duo.cohort.vcf = None
        self.assertIn(MOSAIC_JOINT_CALL_WARNING, node.get_warnings())

    def test_mosaic_needs_no_affected_parent(self):
        self.assertEqual(DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.MOSAIC_PARENT), [])

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
        self.duo.relative_affected = True
        try:
            ids = self._filter_variants(self._make_node(DuoInheritance.ANY_AFFECTED))
            self.assertIn(self.parent_only_v.pk, ids)
        finally:
            self.duo.relative_affected = False

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
        self.duo.relative_affected = True
        try:
            self.assertEqual(DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.DOMINANT), [])
        finally:
            self.duo.relative_affected = False

    def test_xlinked_with_father_raises_error(self):
        errors = DuoNode.get_duo_inheritance_errors(self.duo_father, DuoInheritance.XLINKED_RECESSIVE)
        self.assertGreater(len(errors), 0)
        self.assertIn("mother", errors[0].lower())

    def test_xlinked_with_mother_no_errors(self):
        self.assertEqual(DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.XLINKED_RECESSIVE), [])

    def test_xlinked_affected_mother_raises_error(self):
        self.duo.relative_affected = True
        try:
            errors = DuoNode.get_duo_inheritance_errors(self.duo, DuoInheritance.XLINKED_RECESSIVE)
            self.assertGreater(len(errors), 0)
            self.assertIn("unaffected mother", errors[0].lower())
        finally:
            self.duo.relative_affected = False

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
        self.assertIn('XLR:', data[DuoInheritance.ALL_RECESSIVE][f'relative_{DuoRelationship.MOTHER}'])

    def test_zygosity_table_all_recessive_father_has_no_xlinked_branch(self):
        data = DuoNode.get_zygosity_table_data()
        self.assertNotIn('XLR:', data[DuoInheritance.ALL_RECESSIVE][f'relative_{DuoRelationship.FATHER}'])

    def test_zygosity_table_all_recessive_sibling_keeps_the_xlinked_branch(self):
        data = DuoNode.get_zygosity_table_data()
        self.assertIn('XLR:', data[DuoInheritance.ALL_RECESSIVE][f'relative_{DuoRelationship.SIBLING}_affected'])

    def test_zygosity_table_xlinked_has_chr_x_other_filter(self):
        data = DuoNode.get_zygosity_table_data()
        self.assertEqual(data[DuoInheritance.XLINKED_RECESSIVE]['other_filters_relative'], "Chr X only")

    def test_zygosity_table_compound_het_says_one_hit_from_the_parent(self):
        data = DuoNode.get_zygosity_table_data()
        entry = data[DuoInheritance.COMPOUND_HET]
        self.assertIn("gene", entry[f'other_filters_relative_{DuoRelationship.MOTHER}'])
        self.assertIn("unphased", entry[f'other_filters_relative_{DuoRelationship.SIBLING}_affected'])

    def test_zygosity_table_recessive_sibling_row_follows_their_affected_status(self):
        """ An affected sibling is homozygous like the proband; an unaffected one just isn't """
        entry = DuoNode.get_zygosity_table_data()[DuoInheritance.RECESSIVE]
        affected = entry[f'relative_{DuoRelationship.SIBLING}_affected']
        unaffected = entry[f'relative_{DuoRelationship.SIBLING}_unaffected']
        self.assertEqual(entry['proband'], affected)
        self.assertNotIn("HOM_ALT", unaffected)

    def test_zygosity_table_any_affected_relative_only_counted_when_affected(self):
        data = DuoNode.get_zygosity_table_data()
        entry = data[DuoInheritance.ANY_AFFECTED]
        self.assertTrue(entry['relative_affected'])
        self.assertNotEqual(entry['relative_affected'], entry['relative_unaffected'])

    def test_zygosity_table_mosaic_thresholds_are_on_the_relative_row_only(self):
        """ The relative's row carries the threshold template the editor fills in; the proband row is
            about the proband, so it says nothing about parental read support """
        entry = DuoNode.get_zygosity_table_data()[DuoInheritance.MOSAIC_PARENT]
        self.assertIn("{alt_reads}", entry["other_filters_relative"])
        self.assertIn("{af}", entry["other_filters_relative"])
        self.assertNotIn("other_filters_proband", entry)

    # ── Sibling duo (#1861) ───────────────────────────────────────────────────

    def test_recessive_affected_sibling_wants_both_homozygous(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.RECESSIVE, duo=self.duo_sibling_affected))
        self.assertIn(self.both_hom_alt_v.pk, ids)
        self.assertNotIn(self.recessive_v.pk, ids)  # a HET sibling doesn't share the phenotype

    def test_recessive_unaffected_sibling_excludes_a_homozygous_sibling(self):
        ids = self._filter_variants(self._make_node(DuoInheritance.RECESSIVE, duo=self.duo_sibling))
        self.assertIn(self.recessive_v.pk, ids)  # a carrier sibling is fine
        self.assertNotIn(self.both_hom_alt_v.pk, ids)

    def test_dominant_unaffected_sibling_is_a_discordant_pair_filter(self):
        """ Proband has it, sibling doesn't - a valid filter, so no 'needs an affected parent' error """
        self.assertEqual(DuoNode.get_duo_inheritance_errors(self.duo_sibling, DuoInheritance.DOMINANT), [])
        ids = self._filter_variants(self._make_node(DuoInheritance.DOMINANT, duo=self.duo_sibling))
        self.assertIn(self.absent_v.pk, ids)
        self.assertNotIn(self.dominant_v.pk, ids)

    def test_xlinked_with_a_sibling_does_not_need_the_mother(self):
        self.assertEqual(
            DuoNode.get_duo_inheritance_errors(self.duo_sibling, DuoInheritance.XLINKED_RECESSIVE), [])

    def test_xlinked_affected_sibling_matches_shared_x_homozygote(self):
        node = self._make_node(DuoInheritance.XLINKED_RECESSIVE, duo=self.duo_sibling_affected)
        ids = self._filter_variants(node)
        self.assertIn(self.x_both_hom_alt_v.pk, ids)
        self.assertNotIn(self.both_hom_alt_v.pk, ids)

    def test_all_recessive_with_a_sibling_keeps_the_xlinked_branch(self):
        node = self._make_node(DuoInheritance.ALL_RECESSIVE, duo=self.duo_sibling_affected)
        self.assertIn("XLR", node._get_method_summary())
        ids = self._filter_variants(node)
        self.assertIn(self.both_hom_alt_v.pk, ids)
        self.assertIn(self.x_both_hom_alt_v.pk, ids)

    def test_compound_het_affected_sibling_is_het_on_both_hits(self):
        handler = self._make_node(DuoInheritance.COMPOUND_HET,
                                  duo=self.duo_sibling_affected)._inheritance_factory()
        self.assertEqual(handler._mum_but_not_dad(), handler._dad_but_not_mum())
        self.assertEqual({Zygosity.HET}, handler._mum_but_not_dad()[0])

    def test_compound_het_unaffected_sibling_is_unconstrained(self):
        handler = self._make_node(DuoInheritance.COMPOUND_HET, duo=self.duo_sibling)._inheritance_factory()
        self.assertEqual(set(), handler._mum_but_not_dad()[0])

    def test_compound_het_with_a_sibling_warns_it_is_unphased(self):
        node = self._make_node(DuoInheritance.COMPOUND_HET, duo=self.duo_sibling_affected)
        self.assertIn(COMP_HET_SIBLING_UNPHASED, node.get_warnings())

    def test_compound_het_with_a_parent_does_not_warn(self):
        self.assertEqual(self._make_node(DuoInheritance.COMPOUND_HET).get_warnings(), [])

    def test_parent_only_modes_error_on_a_sibling_duo(self):
        for inheritance in (DuoInheritance.ABSENT_IN_PARENT, DuoInheritance.MOSAIC_PARENT):
            errors = DuoNode.get_duo_inheritance_errors(self.duo_sibling, inheritance)
            self.assertEqual(1, len(errors), f"{inheritance} should need a parent")
            self.assertIn("needs a parent", errors[0])
