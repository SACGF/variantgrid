"""
    Sample.delete_internal_data wipes a sample's own stats rows so a VCF can be
    reloaded in place. The cohort-level (sample IS NULL) rows belong to the CGC
    and must survive.
"""
from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.models import (
    CohortGenotypeClinVarAnnotationStats,
    CohortGenotypeGeneAnnotationStats,
    CohortGenotypeVariantAnnotationStats,
)
from snpdb.models import GenomeBuild, SampleLocusCount, SampleStatsCodeVersion
from snpdb.models.models_cohort_stats import CohortGenotypeStats
from snpdb.models.models_enums import ImportStatus
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort


class SampleDeleteInternalDataTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="delete_internal_data_user")[0]
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.cohort = create_fake_cohort(cls.user, cls.genome_build)
        cls.cgc = cls.cohort.cohort_genotype_collection
        cls.sample = cls.cohort.get_samples().first()
        cls.code_version = SampleStatsCodeVersion.objects.create(
            name="test", version=99, code_git_hash="deleteinternaldata")

    def _create_stats(self, sample):
        common = {
            "cohort_genotype_collection": self.cgc,
            "sample": sample,
            "filter_key": None,
            "passing_filter": False,
            "code_version": self.code_version,
        }
        CohortGenotypeStats.objects.create(import_status=ImportStatus.SUCCESS, **common)
        CohortGenotypeVariantAnnotationStats.objects.create(
            variant_annotation_version=self.annotation_version.variant_annotation_version, **common)
        CohortGenotypeGeneAnnotationStats.objects.create(
            gene_annotation_version=self.annotation_version.gene_annotation_version, **common)
        CohortGenotypeClinVarAnnotationStats.objects.create(
            clinvar_version=self.annotation_version.clinvar_version, **common)

    def test_deletes_per_sample_stats_and_keeps_aggregates(self):
        self._create_stats(self.sample)
        self._create_stats(None)
        SampleLocusCount.objects.create(sample=self.sample, locus_count=1, count=5)

        self.sample.delete_internal_data()

        for klass in (CohortGenotypeStats, CohortGenotypeVariantAnnotationStats,
                      CohortGenotypeGeneAnnotationStats, CohortGenotypeClinVarAnnotationStats):
            self.assertFalse(klass.objects.filter(sample=self.sample).exists(),
                             f"{klass.__name__} per-sample row deleted")
            self.assertTrue(klass.objects.filter(cohort_genotype_collection=self.cgc,
                                                 sample__isnull=True).exists(),
                            f"{klass.__name__} aggregate row kept")
        self.assertFalse(SampleLocusCount.objects.filter(sample=self.sample).exists())
