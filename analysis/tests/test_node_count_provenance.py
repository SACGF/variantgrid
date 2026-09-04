"""
Node count provenance - a node declares the mutable tables its query reads, and the versions of those
sources are recorded with the count. Consumers use that to tell an exact count from an advisory one.
"""
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import Analysis
from analysis.models.nodes.analysis_node import AnalysisNode, NodeVersion
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.models.nodes.sources.sample_node import SampleNode
from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.django_partition import temporary_db_table
from snpdb.models import BuiltInFilters, GenomeBuild, VariantZygosityCountCollection
from snpdb.models.models_cohort import CohortGenotype, CohortGenotypeCollection
from snpdb.models.models_enums import CohortGenotypeCollectionType
from snpdb.tests.utils.fake_cohort_data import create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestNodeCountProvenance(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='node_count_provenance_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.vzcc = VariantZygosityCountCollection.objects.get_or_create(
            name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)[0]

        cls.trio = create_fake_trio(cls.user, cls.grch37)
        cls.cohort = cls.trio.cohort
        cls.proband = cls.cohort.cohortsample_set.get(sample__name='proband').sample
        cls.cgc = CohortGenotypeCollection.objects.get(
            cohort=cls.cohort, cohort_version=cls.cohort.version,
            collection_type=CohortGenotypeCollectionType.UNCOMMON)

        # 3 variants the proband is called on, so a SampleNode has a small non-zero count
        cls.variants = []
        for position in (1000, 2000, 3000):
            variant = slowly_create_test_variant("1", position, "A", "T", cls.grch37)
            cls._add_genotype(variant, "E..")
            cls.variants.append(variant)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)

    @classmethod
    def _add_genotype(cls, variant, samples_zygosity):
        n = len(samples_zygosity)
        with temporary_db_table(CohortGenotype, cls.cgc.get_partition_table()):
            CohortGenotype.objects.create(
                collection=cls.cgc, variant=variant,
                ref_count=samples_zygosity.count('R'),
                het_count=samples_zygosity.count('E'),
                hom_count=samples_zygosity.count('O'),
                unk_count=samples_zygosity.count('U'),
                samples_zygosity=samples_zygosity,
                samples_allele_depth=[20] + [0] * (n - 1),
                samples_allele_frequency=[0.5] + [0.0] * (n - 1),
                samples_read_depth=[40] + [0] * (n - 1),
                samples_genotype_quality=[30] * n,
                samples_phred_likelihood=[0] * n,
            )

    def _load(self, node):
        """ Run the node's counts as the load pipeline would """
        status, count = node.node_counts()
        node.update(status=status, count=count)
        node.status = status
        node.count = count
        return node

    def _source_node(self):
        return self._load(SampleNode.objects.create(analysis=self.analysis, sample=self.proband))

    def _population_child(self, use_internal_counts):
        pop = PopulationNode.objects.create(analysis=self.analysis, use_internal_counts=use_internal_counts,
                                            max_samples=100)
        pop.add_parent(self._source_node())
        pop._cached_parents = None
        pop.save()
        return PopulationNode.objects.get(pk=pop.pk)

    # ── Declaring live data sources ───────────────────────────────────────────

    def test_population_node_internal_counts_declares_zygosity_collection(self):
        pop = self._population_child(use_internal_counts=True)
        self.assertEqual({self.vzcc.live_source_key: self.vzcc.data_version},
                         pop._get_live_data_sources())

    def test_population_node_without_internal_counts_is_deterministic(self):
        pop = self._population_child(use_internal_counts=False)
        self.assertEqual({}, pop._get_live_data_sources())
        self.assertEqual({}, pop.get_live_data_sources(), "Deterministic parent contributes nothing")

    def test_child_inherits_parent_live_data_sources(self):
        pop = self._load(self._population_child(use_internal_counts=True))
        child = PopulationNode.objects.create(analysis=self.analysis, use_internal_counts=False)
        child.add_parent(pop)
        child._cached_parents = None
        child.save()
        child = PopulationNode.objects.get(pk=child.pk)

        self.assertEqual({}, child._get_live_data_sources(), "Child has no live source of its own")
        self.assertIn(self.vzcc.live_source_key, child.get_live_data_sources())

    def test_load_records_live_data_sources_on_node_version(self):
        pop = self._load(self._population_child(use_internal_counts=True))
        node_version = NodeVersion.objects.get(node=pop, version=pop.version)
        self.assertEqual({self.vzcc.live_source_key: self.vzcc.data_version},
                         node_version.live_data_sources)
        self.assertFalse(PopulationNode.objects.get(pk=pop.pk).count_is_deterministic)

    # ── Storing the exact PK set for small nodes ──────────────────────────────

    def test_small_node_stores_its_variant_ids(self):
        node = self._source_node()
        node_version = NodeVersion.objects.get(node=node, version=node.version)

        expected_pks = set(node.get_queryset().values_list("pk", flat=True))
        self.assertEqual(expected_pks, set(node_version.variant_ids))
        self.assertEqual(node_version.counts[BuiltInFilters.TOTAL], len(node_version.variant_ids))
        self.assertEqual(node.count, len(node_version.variant_ids))

    @override_settings(ANALYSIS_NODE_STORE_ID_SIZE_MAX=0)
    def test_large_node_stores_no_variant_ids(self):
        node = self._source_node()
        self.assertIsNone(NodeVersion.objects.get(node=node, version=node.version).variant_ids)
        self.assertIsNone(AnalysisNode.get_cached_node_pks(node))
        self.assertIsNone(AnalysisNode.get_small_parent_arg_q_dict(node),
                          "Caller falls back to the parent subquery")

    def test_small_node_substitutes_explicit_pks(self):
        node = self._source_node()
        arg_q_dict = AnalysisNode.get_small_parent_arg_q_dict(node)
        (q,) = arg_q_dict[None].values()
        self.assertEqual(set(node.get_queryset().values_list("pk", flat=True)),
                         set(q.children[0][1]))

    # ── Provenance-scoped mismatch checks ─────────────────────────────────────

    def test_count_mismatch_raises_on_deterministic_node(self):
        node = self._source_node()
        with self.assertRaises(ValueError):
            node._raise_or_warn_count_mismatch("41 != 42")

    def test_count_mismatch_warns_on_live_source_node(self):
        pop = self._load(self._population_child(use_internal_counts=True))
        with self.assertLogs(level="WARNING") as logs:
            pop._raise_or_warn_count_mismatch("41 != 42")
        self.assertIn(self.vzcc.live_source_key, "\n".join(logs.output))
