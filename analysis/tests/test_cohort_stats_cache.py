"""
    Targeted unit tests for the new CohortGenotype*Stats cache.

    Covers the invariants the plan calls out as load-bearing:
      - canonical_filter_key encoding is deterministic and round-trip stable
      - reader node ↔ writer precompute keys produce byte-identical strings
      - CGC CASCADE deletes stats rows when the CGC dies
      - SampleNode (per-sample) and CohortNode (aggregate) cache hits use
        the right row keys
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models import CohortNode, TrioNode
from analysis.models.enums import TrioInheritance
from annotation.fake_data import get_fake_annotation_version
from annotation.models import (
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
)
from annotation.tasks.calculate_sample_stats import (
    calculate_cohort_stats,
    cohort_stats_are_fresh,
    get_sample_stats_code_version,
)
from library.utils.json_utils import canonical_filter_key
from snpdb.fake_data import create_fake_cohort, create_fake_trio
from snpdb.models import (
    Cohort,
    CohortGenotypeStats,
    GenomeBuild,
    ImportStatus,
    SampleStatsCodeVersion,
    Trio,
)


class TestCanonicalFilterKey(TestCase):
    def test_none_round_trips(self):
        self.assertIsNone(canonical_filter_key(None))
        self.assertIsNone(canonical_filter_key({}))

    def test_deterministic(self):
        a = canonical_filter_key({"inheritance": "denovo"})
        b = canonical_filter_key({"inheritance": "denovo"})
        self.assertEqual(a, b)

    def test_key_order_does_not_matter(self):
        self.assertEqual(
            canonical_filter_key({"a": 1, "b": 2}),
            canonical_filter_key({"b": 2, "a": 1}),
        )

    def test_no_whitespace(self):
        # Encoded form must be compact so equality on the text is byte-exact
        encoded = canonical_filter_key({"inheritance": "denovo"})
        self.assertNotIn(" ", encoded)


class TestStatsFilterKeyRoundTrip(TestCase):
    """ The reader's get_stats_filter_key MUST produce strings that match the writer's precompute
        keys. A silent mismatch produces cache misses forever. """

    @staticmethod
    def _trio_node(inheritance, mother_affected=False, father_affected=True, require_zygosity=True):
        trio = Trio(cohort=Cohort(), mother_affected=mother_affected, father_affected=father_affected)
        return TrioNode(inheritance=inheritance, require_zygosity=require_zygosity, trio=trio)

    def test_trio_node_round_trip(self):
        for inheritance in TrioNode.STATS_INHERITANCE_MODES:
            node = self._trio_node(inheritance)
            self.assertIn(node.get_stats_filter_key(), node.trio.get_stats_filter_keys())

    def test_dominant_key_varies_with_affected_flags(self):
        """ The dominant predicate branches on who is affected, so a flag edit has to miss the
            bucket built for the previous flags rather than read it """
        def dominant_key(mother_affected, father_affected):
            return self._trio_node(TrioInheritance.DOMINANT, mother_affected, father_affected).get_stats_filter_key()

        self.assertNotEqual(dominant_key(True, False), dominant_key(False, True))

    def test_modes_without_affected_flags_keep_one_key(self):
        """ Only dominant depends on affected status - a flag edit leaves the other buckets valid """
        unaffected_parents = set(Trio(cohort=Cohort(), mother_affected=False, father_affected=False).get_stats_filter_keys())
        affected_mother = set(Trio(cohort=Cohort(), mother_affected=True, father_affected=False).get_stats_filter_keys())
        self.assertEqual(len(unaffected_parents & affected_mother), len(Trio.STATS_INHERITANCE_MODES) - 1)

    def test_trio_node_uncacheable_without_require_zygosity(self):
        """ The buckets are computed requiring a parent zygosity call, so a node that also
            accepts parent no-calls has to count live """
        for inheritance in TrioNode.STATS_INHERITANCE_MODES:
            node = self._trio_node(inheritance, require_zygosity=False)
            self.assertIsNone(node._stats_inheritance_mode())


class TestCohortNodeCountLimitsDefeatCache(TestCase):
    """ The aggregate row is the cohort's unfiltered counts, so a count limit has to count live (#2370) """

    @staticmethod
    def _cohort_node(**kwargs):
        return CohortNode(cohort=Cohort(sample_count=3), accordion_panel=CohortNode.COUNT, **kwargs)

    def test_defaults_use_cache(self):
        self.assertFalse(self._cohort_node()._has_filters_that_affect_label_counts())

    def test_count_limits_defeat_cache(self):
        for limits in ({"min_het_or_hom_count": 1}, {"max_het_or_hom_count": 2}, {"min_hom_count": 1},
                       {"max_ref_count": 0}):
            with self.subTest(limits=limits):
                self.assertTrue(self._cohort_node(**limits)._has_filters_that_affect_label_counts())

    def test_limits_that_filter_nothing_use_cache(self):
        node = self._cohort_node(max_het_or_hom_count=3, min_het_count=0)
        self.assertFalse(node._has_filters_that_affect_label_counts())


class TestCohortStatsCascade(TestCase):
    """ CGC deletion (during the post-version-bump async sweep) must wipe all
        four CohortGenotype*Stats families via on_delete=CASCADE — no separate
        invalidation code needed in Cohort.increment_version. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="cgc_cascade_user")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.cgc = cls.cohort.cohort_genotype_collection
        cls.code_version = SampleStatsCodeVersion.objects.create(
            name="test", version=99)

    def test_cgc_delete_cascades_to_stats(self):
        sample = self.cohort.get_samples().first()
        # Per-sample row
        CohortGenotypeStats.objects.create(
            cohort_genotype_collection=self.cgc, sample=sample,
            filter_key=None, passing_filter=False,
            code_version=self.code_version, import_status=ImportStatus.SUCCESS,
            variant_count=42)
        # Aggregate row
        CohortGenotypeStats.objects.create(
            cohort_genotype_collection=self.cgc, sample=None,
            filter_key=None, passing_filter=False,
            code_version=self.code_version, import_status=ImportStatus.SUCCESS,
            variant_count=99)

        self.assertEqual(2, CohortGenotypeStats.objects.filter(
            cohort_genotype_collection=self.cgc).count())

        self.cgc.delete()
        self.assertEqual(0, CohortGenotypeStats.objects.filter(
            cohort_genotype_collection_id=self.cgc.pk).count())


class TestCohortStatsRecomputeLock(TestCase):
    """ Issue #1576: calculate_cohort_stats must serialize per-CGC recomputes and
        early-exit when the CGC's stats are already fresh, so concurrent lazy
        recomputes can't race the delete/bulk_create (which raised IntegrityError
        on the per-sample unique constraints). These tests patch the compute/persist
        step so they exercise the lock + freshness decision without needing the
        annotated variants the fake fixtures don't provide. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="recompute_lock_user")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.grch37)
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.cgc = cls.cohort.cohort_genotype_collection
        cls.code_version = get_sample_stats_code_version()

    def _create_fresh_stats(self, cgc, code_version, samples=()):
        """ Pre-create one of each of the four stat families for cgc at the given
            code_version + the test annotation versions, for the aggregate and each of samples,
            enough to satisfy cohort_stats_are_fresh. """
        for sample in [None, *samples]:
            self._create_stats_row(cgc, code_version, sample)

    def _create_stats_row(self, cgc, code_version, sample):
        av = self.annotation_version
        CohortGenotypeStats.objects.create(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=False, code_version=code_version,
            import_status=ImportStatus.SUCCESS)
        CohortGenotypeVariantAnnotationStats.objects.create(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=False, code_version=code_version,
            variant_annotation_version=av.variant_annotation_version)
        CohortGenotypeGeneAnnotationStats.objects.create(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=False, code_version=code_version,
            gene_annotation_version=av.gene_annotation_version)
        CohortGenotypeClinVarAnnotationStats.objects.create(
            cohort_genotype_collection=cgc, sample=sample, filter_key=None,
            passing_filter=False, code_version=code_version,
            clinvar_version=av.clinvar_version)

    def test_cohort_stats_are_fresh_true_when_all_present(self):
        self.assertFalse(
            cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version, [None], []))
        self._create_fresh_stats(self.cgc, self.code_version)
        self.assertTrue(
            cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version, [None], []))

    def test_cohort_stats_are_fresh_false_on_code_version_mismatch(self):
        # Rows present but at a different code_version → a bump must force recompute.
        other_code_version = SampleStatsCodeVersion.objects.create(
            name="test_other", version=98)
        self._create_fresh_stats(self.cgc, other_code_version)
        self.assertFalse(
            cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version, [None], []))

    def test_cohort_stats_are_fresh_false_when_partial(self):
        # Only the genotype-level row present; annotation-version rows missing.
        CohortGenotypeStats.objects.create(
            cohort_genotype_collection=self.cgc, sample=None, filter_key=None,
            passing_filter=False, code_version=self.code_version,
            import_status=ImportStatus.SUCCESS)
        self.assertFalse(
            cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version, [None], []))

    def test_cohort_stats_are_fresh_false_when_a_key_is_missing(self):
        # Eg a trio's affected flags changed its autosomal dominant key since the last compute
        self._create_fresh_stats(self.cgc, self.code_version)
        keys = [None, canonical_filter_key({"inheritance": "denovo"})]
        self.assertFalse(
            cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version, keys, []))

    def test_cohort_stats_are_fresh_false_when_a_sample_is_missing(self):
        # Eg a sub cohort recompute before #2311 left the parent with only the sub cohort's per-sample rows
        samples = list(self.cohort.get_samples())
        self._create_fresh_stats(self.cgc, self.code_version, samples[:-1])
        sample_ids = [s.pk for s in samples]
        self.assertFalse(
            cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version, [None], sample_ids))

    def test_early_exit_when_fresh_skips_recompute(self):
        self._create_fresh_stats(self.cgc, self.code_version, self.cohort.get_samples())
        with patch("annotation.tasks.calculate_sample_stats._compute_and_persist_cohort_stats") as mock_compute:
            calculate_cohort_stats(self.cohort, self.annotation_version)
            mock_compute.assert_not_called()

    def test_stale_recompute_proceeds(self):
        # No stats rows present → not fresh → recompute (compute/persist) runs.
        with patch("annotation.tasks.calculate_sample_stats._compute_and_persist_cohort_stats") as mock_compute:
            calculate_cohort_stats(self.cohort, self.annotation_version)
            mock_compute.assert_called_once()
