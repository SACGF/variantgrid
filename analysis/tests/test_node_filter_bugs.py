"""
Tests for analysis node filters where a count of 0 is a real user choice rather
than "unset" - these all used to be lost to truthiness checks.
"""

from django.conf import settings
from django.test import TestCase, override_settings

from analysis.enums import GroupOperation
from analysis.models import AllVariantsNode, Analysis
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.models.nodes.sources.sample_node import SampleNode
from analysis.tests.utils import AnalysisSetupMixin, create_annotated_variants, node_keeps
from snpdb.fake_data import create_fake_cohort
from snpdb.models.models_cohort import CohortGenotype
from snpdb.models.models_zygosity_counts import VariantZygosityCount, VariantZygosityCountCollection
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

# ---------------------------------------------------------------------------
# AbstractZygosityCountNode - max_count=0
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestZygosityCountZero(AnalysisSetupMixin, TestCase):
    """ max_count=0 means "seen in no samples" - a real filter, not an unset one. """

    def _node_with_max_zero(self):
        node = AllVariantsNode(analysis=self.analysis)
        node.max_het_or_hom_count = 0
        return node

    def test_max_count_zero_produces_lte_filter(self):
        """ The counts come from the global germline collection """
        vzcc = VariantZygosityCountCollection.objects.get_or_create(
            name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)[0]
        variants = {}
        for i, (name, het_count) in enumerate([("seen_nowhere", 0), ("seen_once", 1)]):
            variant = slowly_create_test_variant("1", 1000 + i * 100, "A", "T", self.grch37)
            VariantZygosityCount.objects.create(variant=variant, collection=vzcc, ref_count=2, het_count=het_count,
                                                hom_count=0)
            variants[variant.pk] = name

        node = self._node_with_max_zero()
        node.save()
        filtered_pks = node.get_queryset().filter(pk__in=variants).values_list("pk", flat=True)
        self.assertEqual({"seen_nowhere"}, {variants[pk] for pk in filtered_pks})

    def test_description_shows_max_zero(self):
        """_get_zygosity_count_description() must include '0' when max_count=0."""
        node = self._node_with_max_zero()
        desc = node._get_zygosity_count_description()
        self.assertIn("0", desc, f"max_count=0 disappeared from description. Got: {desc!r}")


# ---------------------------------------------------------------------------
# PopulationNode.modifies_parents() - gnomad_hom_alt_max=0
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPopulationNodeModifiesParents(AnalysisSetupMixin, TestCase):
    """
    gnomad_hom_alt_max=0 (remove variants with ANY gnomAD hom-alt carrier) is a real
    filter, not "unset" - modifies_parents() has to report it or _get_node_q()'s Q is
    never applied to the queryset.
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vav = cls.analysis.annotation_version.variant_annotation_version
        cls.variants = create_annotated_variants(vav, {
            name: {"gnomad_hom_alt": count, "gnomad_hemi_count": count}
            for name, count in [("none", 0), ("one", 1), ("not_in_gnomad", None)]
        })

    def _node(self, **kwargs):
        """Return an unsaved PopulationNode with 100% (no population AF filter) by default."""
        defaults = {"percent": PopulationNode.EVERYTHING, "show_gnomad_filtered": True}
        defaults.update(kwargs)
        return PopulationNode(analysis=self.analysis, **defaults)

    def _filtered(self, node) -> set[str]:
        return node_keeps(node, self.variants)

    def test_gnomad_hom_alt_max_zero_modifies_parents(self):
        node = self._node(gnomad_hom_alt_max=0)
        self.assertTrue(node.modifies_parents())

    def test_gnomad_hom_alt_max_zero_removes_hom_alt_carriers(self):
        node = self._node(gnomad_hom_alt_max=0)
        self.assertEqual({"none", "not_in_gnomad"}, self._filtered(node))

    def test_gnomad_hom_alt_max_none_does_not_modify_alone(self):
        node = self._node(gnomad_hom_alt_max=None)
        self.assertFalse(node.modifies_parents())

    def test_gnomad_hom_alt_max_positive_modifies_parents(self):
        node = self._node(gnomad_hom_alt_max=5)
        self.assertTrue(node.modifies_parents())

    def test_gnomad_hemi_count_max_zero_filters_hemi_count(self):
        node = self._node(gnomad_hemi_count_max=0)
        self.assertTrue(node.modifies_parents())
        self.assertEqual({"none", "not_in_gnomad"}, self._filtered(node))


# ---------------------------------------------------------------------------
# PopulationNode group operation Q structure
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestPopulationNodeGroupOperationQ(AnalysisSetupMixin, TestCase):
    """
    PopulationNode uses *inverted* operators for ANY/ALL (see comment in source):

      ANY (strict)  — remove if above cutoff in ANY database
      ALL (lenient) — remove only if above cutoff in ALL databases

    These tests pin which variants each keeps, so the inversion isn't "corrected" away.
    """
    RARE = 0.001
    COMMON = 0.05

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        vav = cls.analysis.annotation_version.variant_annotation_version
        cls.variants = create_annotated_variants(vav, {
            name: {"gnomad_af": gnomad_af, "af_1kg": other_af, "af_uk10k": other_af}
            for name, gnomad_af, other_af in [("rare_everywhere", cls.RARE, cls.RARE),
                                              ("common_in_gnomad_only", cls.COMMON, cls.RARE),
                                              ("common_everywhere", cls.COMMON, cls.COMMON)]
        })

    def _filtered(self, group_operation) -> set[str]:
        """ gnomad_af + af_1kg + af_uk10k at 1%. Saved so populationnodegnomadpopulation_set.all() can run """
        node = PopulationNode.objects.create(
            analysis=self.analysis,
            percent=1.0,
            group_operation=group_operation,
            gnomad_af=True,
            gnomad_popmax_af=False,
            af_1kg=True,
            af_uk10k=True,
            topmed_af=False,
        )
        return node_keeps(node, self.variants)

    def test_any_group_operation_combines_with_and(self):
        """ GroupOperation.ANY - common in one database is enough to remove it """
        self.assertEqual({"rare_everywhere"}, self._filtered(GroupOperation.ANY))

    def test_all_group_operation_combines_with_or(self):
        """ GroupOperation.ALL - only removed when common in every database """
        self.assertEqual({"rare_everywhere", "common_in_gnomad_only"}, self._filtered(GroupOperation.ALL))


# ---------------------------------------------------------------------------
# AbstractZygosityCountNode - bounds that match every sample
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestZygosityCountNoOpBounds(AnalysisSetupMixin, TestCase):
    """ A min of 0 or a max of every sample filters nothing, so it stays out of the query -
        that's what keeps AllVariantsNode from joining the global zygosity counts (#24) """

    MAX_SAMPLES = 10

    def _node(self, **kwargs):
        node = AllVariantsNode(analysis=self.analysis, **kwargs)
        # Prime the cached_property behind zygosity_count_max_samples
        node.__dict__["num_samples_for_build"] = self.MAX_SAMPLES
        return node

    def test_default_min_zero_produces_no_filter(self):
        """ min_het_or_hom_count defaults to 0 - every variant is >= 0 """
        self.assertEqual({}, self._node().get_zygosity_count_arg_q_dict())

    def test_max_at_sample_count_produces_no_filter(self):
        node = self._node(max_het_or_hom_count=self.MAX_SAMPLES)
        self.assertEqual({}, node.get_zygosity_count_arg_q_dict())

    def test_max_below_sample_count_filters(self):
        node = self._node(max_het_or_hom_count=self.MAX_SAMPLES - 1)
        self.assertTrue(node.get_zygosity_count_arg_q_dict())

    def test_no_op_bounds_left_out_of_description(self):
        node = self._node(min_het_or_hom_count=0, max_het_or_hom_count=self.MAX_SAMPLES,
                          min_het_count=2)
        self.assertEqual("Het >= 2", node._get_zygosity_count_description())


# ---------------------------------------------------------------------------
# SampleNode thresholds vs the columns the VCF carries - max_pl=0
# ---------------------------------------------------------------------------

@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestSampleNodeThresholdColumns(AnalysisSetupMixin, TestCase):
    """ PL<=0 is a real filter, but a caller that reports no PL (TSO 500 fusions) leaves the column
        null on every row, so applying it emptied the node - while the cached stats lookup read
        max_pl=0 as unset and answered with the sample's whole count """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.cohort = create_fake_cohort(cls.analysis.user, cls.grch37)
        cls.sample = cls.cohort.get_samples().first()

    def _node(self, **kwargs):
        return SampleNode(analysis=self.analysis, sample=self.sample, **kwargs)

    def _genotype(self, position: int, proband_pl) -> int:
        """ A het proband call - proband_pl of None is a caller that reports no PL at all """
        variant = slowly_create_test_variant("1", position, "A", "T", self.grch37)
        CohortGenotype.objects.create(
            collection=self.cohort.cohort_genotype_collection, variant=variant, ref_count=2, het_count=1,
            samples_zygosity="ERR", samples_allele_depth=[20] * 3, samples_allele_frequency=[0.5] * 3,
            samples_read_depth=[30] * 3, samples_genotype_quality=[30] * 3,
            samples_phred_likelihood=None if proband_pl is None else [proband_pl, 0, 0])
        return variant.pk

    @staticmethod
    def _filtered_pks(node, pks) -> set[int]:
        node.save()
        return set(node.get_queryset().filter(pk__in=pks).values_list("pk", flat=True))

    def test_max_pl_left_out_when_vcf_has_no_pl(self):
        no_pl_pk = self._genotype(1000, None)
        node = self._node(max_pl=0)
        self.assertEqual({no_pl_pk}, self._filtered_pks(node, [no_pl_pk]))
        self.assertFalse(node._has_filters_that_affect_label_counts())

    def test_max_pl_filters_when_vcf_has_pl(self):
        self.sample.vcf.phred_likelihood_field = "PL"
        self.sample.vcf.save()
        pl_0_pk = self._genotype(1000, 0)
        pl_30_pk = self._genotype(1100, 30)
        node = self._node(max_pl=0)
        self.assertEqual({pl_0_pk}, self._filtered_pks(node, [pl_0_pk, pl_30_pk]))
        self.assertTrue(node._has_filters_that_affect_label_counts())

    def test_min_ad_filters_when_vcf_has_ad(self):
        node = self._node(min_ad=5)
        self.assertTrue(node._has_filters_that_affect_label_counts())

    def test_min_dp_left_out_when_vcf_has_no_dp(self):
        node = self._node(min_dp=5)
        self.assertFalse(node._has_filters_that_affect_label_counts())

    def test_method_summary_only_shows_applied_thresholds(self):
        node = self._node(min_ad=5, min_dp=10, max_pl=0)
        summary = node._get_sample_method_summary(self.sample)
        self.assertIn("AD>=5", summary)
        self.assertNotIn("DP>=", summary)
        self.assertNotIn("PL<=", summary)
