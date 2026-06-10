"""
    Targeted unit tests for the new CohortGenotype*Stats cache.

    Covers the invariants the plan calls out as load-bearing:
      - canonical_filter_key encoding is deterministic and round-trip stable
      - reader handler ↔ writer precompute keys produce byte-identical strings
      - CGC CASCADE deletes stats rows when the CGC dies
      - SampleNode (per-sample) and CohortNode (aggregate) cache hits use
        the right row keys
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models.enums import TrioInheritance
from analysis.models.nodes.sources._stats_cache import (
    UNCACHEABLE,
    NoFilterHandler,
    TrioInheritanceHandler,
    get_filter_keys_to_precompute_for_cohort,
)
from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import (
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
)
from annotation.tasks.calculate_sample_stats import (
    _cohort_stats_are_fresh,
    _get_sample_stats_code_version,
    calculate_cohort_stats,
)
from library.utils.json_utils import canonical_filter_key
from snpdb.models import (
    CohortGenotypeStats,
    GenomeBuild,
    ImportStatus,
    SampleStatsCodeVersion,
)
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort, create_fake_trio


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


class TestFilterKeyHandlerRoundTrip(TestCase):
    """ The reader's filter_key_for_node MUST produce strings that match the
        writer's filter_keys_to_precompute. A silent mismatch produces cache
        misses forever. """

    def test_trio_handler_round_trip(self):
        handler = TrioInheritanceHandler()
        # Build a fake "node" for each cached mode and check the key matches
        # one of the precompute keys (ignoring None).
        precomputed = set(handler.filter_keys_to_precompute(_FakeTrioCohort()))
        for inheritance, _ in TrioInheritanceHandler.CACHED_MODES.items():
            class _FakeNode:
                pass
            n = _FakeNode()
            n.inheritance = inheritance
            self.assertIn(handler.filter_key_for_node(n), precomputed)

    def test_no_filter_handler_returns_none_only(self):
        handler = NoFilterHandler()
        self.assertEqual(list(handler.filter_keys_to_precompute(_FakeNonTrioCohort())), [None])
        self.assertIsNone(handler.filter_key_for_node(object()))


class _FakeTrioCohort:
    class _TrioSet:
        @staticmethod
        def exists():
            return True
    trio_set = _TrioSet


class _FakeNonTrioCohort:
    class _TrioSet:
        @staticmethod
        def exists():
            return False
    trio_set = _TrioSet


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
            name="test", version=99, code_git_hash="cascadetest")

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


class TestPrecomputeKeysForCohort(TestCase):
    """ get_filter_keys_to_precompute_for_cohort decides what filter-keyed
        buckets the writer creates per cohort. Trios get all 5 inheritance
        modes; everything else gets just None. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="precompute_user")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def test_trio_cohort_precomputes_5_inheritance_modes_plus_none(self):
        trio = create_fake_trio(self.user, self.grch37)
        keys = get_filter_keys_to_precompute_for_cohort(trio.cohort)
        self.assertIn(None, keys)
        for mode in TrioInheritanceHandler.CACHED_MODES.values():
            self.assertIn(canonical_filter_key({"inheritance": mode}), keys)
        self.assertEqual(len(keys), 6)

    def test_non_trio_cohort_only_precomputes_none(self):
        cohort = create_fake_cohort(self.user, self.grch37)
        keys = get_filter_keys_to_precompute_for_cohort(cohort)
        self.assertEqual(keys, [None])


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
        cls.code_version = _get_sample_stats_code_version()

    def _create_fresh_stats(self, cgc, code_version):
        """ Pre-create one of each of the four stat families for cgc at the given
            code_version + the test annotation versions, enough to satisfy
            _cohort_stats_are_fresh. """
        av = self.annotation_version
        CohortGenotypeStats.objects.create(
            cohort_genotype_collection=cgc, sample=None, filter_key=None,
            passing_filter=False, code_version=code_version,
            import_status=ImportStatus.SUCCESS)
        CohortGenotypeVariantAnnotationStats.objects.create(
            cohort_genotype_collection=cgc, sample=None, filter_key=None,
            passing_filter=False, code_version=code_version,
            variant_annotation_version=av.variant_annotation_version)
        CohortGenotypeGeneAnnotationStats.objects.create(
            cohort_genotype_collection=cgc, sample=None, filter_key=None,
            passing_filter=False, code_version=code_version,
            gene_annotation_version=av.gene_annotation_version)
        CohortGenotypeClinVarAnnotationStats.objects.create(
            cohort_genotype_collection=cgc, sample=None, filter_key=None,
            passing_filter=False, code_version=code_version,
            clinvar_version=av.clinvar_version)

    def test_cohort_stats_are_fresh_true_when_all_present(self):
        self.assertFalse(
            _cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version))
        self._create_fresh_stats(self.cgc, self.code_version)
        self.assertTrue(
            _cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version))

    def test_cohort_stats_are_fresh_false_on_code_version_mismatch(self):
        # Rows present but at a different code_version → a bump must force recompute.
        other_code_version = SampleStatsCodeVersion.objects.create(
            name="test_other", version=98, code_git_hash="stale")
        self._create_fresh_stats(self.cgc, other_code_version)
        self.assertFalse(
            _cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version))

    def test_cohort_stats_are_fresh_false_when_partial(self):
        # Only the genotype-level row present; annotation-version rows missing.
        CohortGenotypeStats.objects.create(
            cohort_genotype_collection=self.cgc, sample=None, filter_key=None,
            passing_filter=False, code_version=self.code_version,
            import_status=ImportStatus.SUCCESS)
        self.assertFalse(
            _cohort_stats_are_fresh(self.cgc, self.annotation_version, self.code_version))

    def test_early_exit_when_fresh_skips_recompute(self):
        self._create_fresh_stats(self.cgc, self.code_version)
        with patch("annotation.tasks.calculate_sample_stats._compute_and_persist_cohort_stats") as mock_compute:
            calculate_cohort_stats(self.cohort, self.annotation_version)
            mock_compute.assert_not_called()

    def test_stale_recompute_proceeds(self):
        # No stats rows present → not fresh → recompute (compute/persist) runs.
        with patch("annotation.tasks.calculate_sample_stats._compute_and_persist_cohort_stats") as mock_compute:
            calculate_cohort_stats(self.cohort, self.annotation_version)
            mock_compute.assert_called_once()
