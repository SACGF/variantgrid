from unittest import mock

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models import GenomeBuild, VariantZygosityCountCollection
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort
from snpdb.variant_zygosity_count import update_variant_zygosity_count_for_vcf


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
