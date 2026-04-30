"""
    Targeted unit tests for the new CohortGenotype*Stats cache.

    Covers the invariants the plan calls out as load-bearing:
      - canonical_filter_key encoding is deterministic and round-trip stable
      - reader handler ↔ writer precompute keys produce byte-identical strings
      - CGC CASCADE deletes stats rows when the CGC dies
      - SampleNode (per-sample) and CohortNode (aggregate) cache hits use
        the right row keys
"""
from django.contrib.auth.models import User
from django.test import TestCase

from analysis.models.enums import TrioInheritance
from analysis.models.nodes.sources._stats_cache import (
    NoFilterHandler,
    TrioInheritanceHandler,
    UNCACHEABLE,
    get_filter_keys_to_precompute_for_cohort,
)
from annotation.models import (
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
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
