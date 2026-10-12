"""
The node grid's genotype columns read the node's own cohort genotype join (#2074).

A node between the source and a rare PopulationNode widens the join to the common collection, while the
grid on its own would ask for the rare collection only - so the grid must reuse the node's join rather
than add a second one, or rows the node took from the common collection show blank genotypes.
"""
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.grids import VariantGrid
from analysis.models import Analysis
from analysis.models.enums import ZygosityNodeZygosity
from analysis.models.nodes.filters.population_node import PopulationNode
from analysis.models.nodes.filters.zygosity_node import ZygosityNode
from analysis.models.nodes.sources.sample_node import SampleNode
from annotation.fake_data import get_fake_annotation_version
from annotation.models import AnnotationRun
from annotation.models.models import VariantAnnotation, VariantAnnotationVersion
from library.django_utils import FakeRequest
from library.django_utils.django_partition import temporary_db_table
from snpdb.fake_data import create_fake_trio
from snpdb.models import GenomeBuild
from snpdb.models.models_cohort import (
    CohortGenotype,
    CohortGenotypeCollection,
    CohortGenotypeCommonFilterVersion,
)
from snpdb.models.models_enums import CohortGenotypeCollectionType
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False, ANALYSIS_NODE_STORE_ID_SIZE_MAX=0)
class TestNodeGridGenotypeJoin(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_2074_grid_join')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.annotation_run = AnnotationRun.objects.create()

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(cls.user)
        cls.vav = cls.analysis.annotation_version.variant_annotation_version

        cohort = create_fake_trio(cls.user, cls.grch37).cohort
        cls.proband = cohort.cohortsample_set.get(sample__name='proband').sample
        cls.uncommon_cgc = CohortGenotypeCollection.objects.get(
            cohort=cohort, cohort_version=cohort.version, collection_type=CohortGenotypeCollectionType.UNCOMMON)
        cfv = CohortGenotypeCommonFilterVersion.objects.create(
            gnomad_version="2.1.1", gnomad_af_min=0.05, genome_build=cls.grch37)
        common_cgc = CohortGenotypeCollection.objects.create(
            cohort=cohort, cohort_version=cohort.version, num_samples=cls.uncommon_cgc.num_samples,
            collection_type=CohortGenotypeCollectionType.COMMON, common_filter=cfv)
        cls.uncommon_cgc.common_collection = common_cgc
        cls.uncommon_cgc.save()

        cls.v_rare = slowly_create_test_variant("1", 1000, "A", "T", cls.grch37)
        cls._add_genotype(cls.uncommon_cgc, cls.v_rare, "E..")
        cls._add_annotation(cls.v_rare, 0.001)
        # In the common collection, but with no population frequency the rare filter keeps it
        cls.v_common = slowly_create_test_variant("1", 4000, "A", "T", cls.grch37)
        cls._add_genotype(common_cgc, cls.v_common, "E..")
        cls._add_annotation(cls.v_common, None)

    @classmethod
    def _add_genotype(cls, cgc, variant, samples_zygosity):
        n = len(samples_zygosity)
        with temporary_db_table(CohortGenotype, cgc.get_partition_table()):
            CohortGenotype.objects.create(
                collection=cgc, variant=variant,
                het_count=samples_zygosity.count('E'), ref_count=0, hom_count=0, unk_count=0,
                samples_zygosity=samples_zygosity,
                samples_allele_depth=[20] * n, samples_allele_frequency=[0.5] * n,
                samples_read_depth=[40] * n, samples_genotype_quality=[30] * n,
                samples_phred_likelihood=[0] * n)

    @classmethod
    def _add_annotation(cls, variant, af):
        partition_table = cls.vav.get_partition_table(
            base_table_name=VariantAnnotationVersion.REPRESENTATIVE_TRANSCRIPT_ANNOTATION)
        with temporary_db_table(VariantAnnotation, partition_table):
            VariantAnnotation.objects.create(
                version=cls.vav, variant=variant, annotation_run=cls.annotation_run,
                gnomad_af=af, af_1kg=af, af_uk10k=af,
                predictions_num_pathogenic=0, predictions_num_benign=0)

    @staticmethod
    def _ready(node):
        status, count = node.node_counts()
        node.update(status=status, count=count)
        node.status = status
        node.count = count
        return node

    def _add_child(self, parent, child):
        child.add_parent(parent)
        child._cached_parents = None
        child.save()
        return child

    def test_rare_population_below_intermediate_node(self):
        sample_node = self._ready(SampleNode.objects.create(analysis=self.analysis, sample=self.proband))
        zygosity_node = self._ready(self._add_child(sample_node, ZygosityNode.objects.create(
            analysis=self.analysis, sample=self.proband, zygosity=ZygosityNodeZygosity.HET)))
        population_node = self._add_child(zygosity_node, PopulationNode.objects.create(
            analysis=self.analysis, percent=1.0, gnomad_af=True, gnomad_popmax_af=False,
            af_1kg=True, af_uk10k=True, topmed_af=False))

        grid = VariantGrid(FakeRequest(user=self.user), population_node)
        qs = grid.get_initial_queryset()
        self.assertEqual(str(qs.query).count('JOIN "snpdb_cohortgenotype"'), 1)

        zygosity_alias = self.uncommon_cgc.get_packed_column_alias("samples_zygosity")
        zygosity_by_variant = dict(qs.values_list("pk", zygosity_alias))
        self.assertEqual(zygosity_by_variant, {self.v_rare.pk: "E..", self.v_common.pk: "E.."})
