from unittest import mock

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from patients.enums import Zygosity
from snpdb.fake_data import create_fake_cohort
from snpdb.models import (
    CohortGenotype,
    CohortGenotypeCollection,
    GenomeBuild,
    VariantZygosityCount,
    VariantZygosityCountCollection,
    VariantZygosityCountForSample,
)
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.variant_zygosity_count import (
    create_variant_zygosity_counts,
    update_variant_zygosity_count_for_sample,
    update_variant_zygosity_count_for_vcf,
)


class VariantZygosityCollectionDataVersionTest(TestCase):
    """ data_version is what a node records as the version of the live counts it filtered on, so it has
        to move with the counts - never ahead of a count update that then fails """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='vzcc_data_version_test_user')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, cls.grch37)
        cls.vcf = cls.cohort.vcf

    def _collection(self) -> VariantZygosityCountCollection:
        return VariantZygosityCountCollection.objects.get_or_create(
            name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)[0]

    def _data_version(self, collection) -> int:
        return VariantZygosityCountCollection.objects.get(pk=collection.pk).data_version

    def test_count_update_bumps_data_version(self):
        collection = self._collection()
        before = self._data_version(collection)

        with mock.patch("snpdb.variant_zygosity_count.run_sql"):
            update_variant_zygosity_count_for_vcf(collection, self.vcf, '+')

        self.assertEqual(before + 1, self._data_version(collection))

    def test_failed_count_update_leaves_data_version_unchanged(self):
        collection = self._collection()
        before = self._data_version(collection)

        real_bump = VariantZygosityCountCollection.bump_data_version

        def bump_then_fail(self):
            real_bump(self)
            raise RuntimeError("count update failed after the bump")

        with mock.patch("snpdb.variant_zygosity_count.run_sql"), \
                mock.patch.object(VariantZygosityCountCollection, "bump_data_version", bump_then_fail):
            with self.assertRaises(RuntimeError):
                update_variant_zygosity_count_for_vcf(collection, self.vcf, '+')

        self.assertEqual(before, self._data_version(collection),
                         "Bump is inside the count update transaction, so it rolls back with it")


class VariantZygosityCountForVCFTest(TestCase):
    """ Adding a VCF to the global counts must count it once however often '+' is requested (#2441) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='vzc_for_vcf_test_user')[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, grch37, name="vzc_for_vcf")
        cls.vcf = cls.cohort.vcf
        cls.variant = slowly_create_test_variant("3", 7000, "A", "T", grch37)
        cgc = CohortGenotypeCollection.objects.get(cohort=cls.cohort)
        CohortGenotype.objects.create(collection=cgc, variant=cls.variant,
                                      samples_zygosity=Zygosity.HET + Zygosity.HOM_ALT + Zygosity.HET,
                                      ref_count=0, het_count=2, hom_count=1,
                                      samples_allele_depth=[20] * 3, samples_allele_frequency=[100] * 3,
                                      samples_read_depth=[30] * 3, samples_genotype_quality=[30] * 3,
                                      samples_phred_likelihood=[0] * 3)

    def setUp(self):
        self.collection = VariantZygosityCountCollection.objects.get_or_create(
            name=settings.VARIANT_ZYGOSITY_GLOBAL_COLLECTION)[0]
        create_variant_zygosity_counts()

    def _counts(self) -> tuple[int, int]:
        vzc = VariantZygosityCount.objects.get(collection=self.collection, variant=self.variant)
        return vzc.het_count, vzc.hom_count

    def test_add_twice_counts_once(self):
        update_variant_zygosity_count_for_vcf(self.collection, self.vcf, '+')
        update_variant_zygosity_count_for_vcf(self.collection, self.vcf, '+')
        self.assertEqual((2, 1), self._counts())

    def test_add_after_subtract_counts_again(self):
        update_variant_zygosity_count_for_vcf(self.collection, self.vcf, '+')
        update_variant_zygosity_count_for_vcf(self.collection, self.vcf, '-')
        self.assertEqual((0, 0), self._counts())
        update_variant_zygosity_count_for_vcf(self.collection, self.vcf, '+')
        self.assertEqual((2, 1), self._counts())
        update_variant_zygosity_count_for_vcf(self.collection, self.vcf, '-')
        self.assertEqual((0, 0), self._counts(), "Re-added count can be subtracted again")

    def test_subtract_one_sample_splits_vcf_count(self):
        update_variant_zygosity_count_for_vcf(self.collection, self.vcf, '+')
        father = self.vcf.sample_set.get(name="father")
        update_variant_zygosity_count_for_sample(self.collection, father, '-')
        self.assertEqual((1, 1), self._counts())
        self.assertEqual(3, VariantZygosityCountForSample.objects.filter(collection=self.collection).count())
