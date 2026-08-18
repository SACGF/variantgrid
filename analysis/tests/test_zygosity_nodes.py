from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from analysis.models import AllVariantsNode, Analysis, AnalysisNode, CohortNode
from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import (
    GenomeBuild,
    VariantZygosityCount,
    VariantZygosityCountCollection,
)
from snpdb.models.models_cohort import CohortGenotype, CohortGenotypeCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class TestZygosityNodes(TestCase):
    """
    Trio layout (from create_fake_trio): packed_field_index 0=proband, 1=mother, 2=father.
    The sub cohort is [proband, mother], so only positions 0 + 1 count for it.
    samples_zygosity encoding: E=HET, R=HOM_REF, O=HOM_ALT
    """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()

        user = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version_grch37 = get_fake_annotation_version(cls.grch37)
        cls.trio = create_fake_trio(user, cls.grch37)
        cls.cgc = CohortGenotypeCollection.objects.get(cohort=cls.trio.cohort)
        cls.vzcc = VariantZygosityCountCollection.objects.get_or_create(
            name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)[0]

        samples_list = list(cls.trio.cohort.get_samples()[:2])  # Leave out last sample to be sub-cohort
        cls.sub_cohort = cls.trio.cohort.create_sub_cohort(user, samples_list)

        cls.analysis = Analysis(genome_build=cls.grch37)
        cls.analysis.set_defaults_and_save(user)

        # One of each zygosity across the trio, but only ref+het in the sub cohort
        cls.all_types_v = cls._make_variant(1000, "REO")
        # Also one of each across the trio, and het+hom in the sub cohort
        cls.het_hom_v = cls._make_variant(2000, "EOR")
        cls.all_ref_v = cls._make_variant(3000, "RRR")
        cls.all_het_v = cls._make_variant(4000, "EEE")

    @classmethod
    def _make_variant(cls, position: int, samples_zygosity: str):
        variant = slowly_create_test_variant("3", position, "A", "T", cls.grch37)
        n = len(samples_zygosity)
        counts = {
            "ref_count": samples_zygosity.count('R'),
            "het_count": samples_zygosity.count('E'),
            "hom_count": samples_zygosity.count('O'),
        }
        CohortGenotype.objects.create(
            collection=cls.cgc,
            variant=variant,
            samples_zygosity=samples_zygosity,
            samples_allele_depth=[20] * n,
            samples_allele_frequency=[100] * n,
            samples_read_depth=[30] * n,
            samples_genotype_quality=[30] * n,
            samples_phred_likelihood=[0] * n,
            **counts,
        )
        # AllVariantsNode counts come from the global collection, not the cohort
        VariantZygosityCount.objects.create(variant=variant, collection=cls.vzcc, **counts)
        return variant

    # At least one sample of every zygosity - the all-ref and all-het variants drop out
    ONE_OF_EACH_ZYGOSITY = {
        "min_het_or_hom_count": 1, "max_het_or_hom_count": 99,
        "min_ref_count": 1, "max_ref_count": 99,
        "min_het_count": 1, "max_het_count": 99,
        "min_hom_count": 1, "max_hom_count": 99,
    }

    def test_cohort_node_zyg_filters(self):
        cohort_node = CohortNode.objects.create(analysis=self.analysis, cohort=self.trio.cohort,
                                                accordion_panel=CohortNode.COUNT)
        self._assert_zygosity_filters(cohort_node, {self.all_types_v, self.het_hom_v},
                                      **self.ONE_OF_EACH_ZYGOSITY)

    def test_sub_cohort_node_zyg_filters(self):
        sub_cohort_node = CohortNode.objects.create(analysis=self.analysis, cohort=self.sub_cohort,
                                                accordion_panel=CohortNode.COUNT)
        # Counts are built from positions 0 + 1 only, so "REO" no longer has a hom while "EOR" does
        self._assert_zygosity_filters(sub_cohort_node, {self.het_hom_v},
                                      min_het_count=1, max_het_count=99,
                                      min_hom_count=1, max_hom_count=99)

    def test_all_variants_node_zyg_filters(self):
        all_variants_node = AllVariantsNode(analysis=self.analysis)
        self._assert_zygosity_filters(all_variants_node, {self.all_types_v, self.het_hom_v},
                                      **self.ONE_OF_EACH_ZYGOSITY)

    def _assert_zygosity_filters(self, node: AnalysisNode, expected_variants, **zygosity_counts):
        for field, value in zygosity_counts.items():
            setattr(node, field, value)
        node.save()  # Also makes a whole bunch of needed stuff

        variant_pks = set(node.get_queryset().values_list("pk", flat=True))
        self.assertEqual(variant_pks, {v.pk for v in expected_variants})
