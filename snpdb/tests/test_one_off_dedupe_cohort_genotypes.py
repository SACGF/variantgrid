from io import StringIO
from unittest import mock

from django.contrib.auth.models import User
from django.core.management import call_command
from django.test import TestCase

from library.django_utils.database_utils import run_sql
from library.django_utils.django_partition import temporary_db_table
from snpdb.fake_data import create_fake_cohort, make_cohort_genotype
from snpdb.models.models_cohort import CohortGenotype
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_zygosity_counts import VariantZygosityCount, VariantZygosityCountCollection
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant
from snpdb.variant_zygosity_count import (
    create_variant_zygosity_counts,
    update_variant_zygosity_count_for_vcf,
)


class DedupeCohortGenotypesTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        user = User.objects.get_or_create(username='dedupe_cohort_genotypes_user')[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.identical_variant = slowly_create_test_variant("3", 128198980, 'A', 'T', grch37)
        cls.differing_variant = slowly_create_test_variant("3", 128198990, 'A', 'G', grch37)
        cls.collection = VariantZygosityCountCollection.get_global_germline_counts()
        create_variant_zygosity_counts()

        duplicated_cohort = create_fake_cohort(user, grch37, name="dedupe_duplicated")
        cls.cgc = duplicated_cohort.cohort_genotype_collection
        # Like the partitions that hold duplicates: created without the unique (collection, variant) index
        partition_table = cls.cgc.get_partition_table()
        run_sql(f'ALTER TABLE "{partition_table}" DROP CONSTRAINT "{partition_table}_collection_id_variant_id_key"')
        with temporary_db_table(CohortGenotype, partition_table):
            make_cohort_genotype(cls.cgc, cls.identical_variant, "ERR")
            make_cohort_genotype(cls.cgc, cls.identical_variant, "ERR")
            make_cohort_genotype(cls.cgc, cls.differing_variant, "ERR")
            make_cohort_genotype(cls.cgc, cls.differing_variant, "OOE")
        other_cohort = create_fake_cohort(user, grch37, name="dedupe_other")
        make_cohort_genotype(other_cohort.cohort_genotype_collection, cls.identical_variant, "EEO")

        for cohort in [duplicated_cohort, other_cohort]:
            update_variant_zygosity_count_for_vcf(cls.collection, cohort.vcf, '+')
        # A stored count that disagrees with the rows, which the recount replaces
        VariantZygosityCount.objects.filter(collection=cls.collection).update(ref_count=9, het_count=9, hom_count=9)

    def _genotype_ids(self, variant) -> list[int]:
        qs = CohortGenotype.objects.filter(collection=self.cgc, variant=variant).order_by("pk")
        return list(qs.values_list("pk", flat=True))

    def _counts(self, variant) -> tuple[int, int, int]:
        vzc = VariantZygosityCount.objects.get(collection=self.collection, variant=variant)
        return vzc.ref_count, vzc.het_count, vzc.hom_count

    @mock.patch("snpdb.management.commands.one_off_dedupe_cohort_genotypes.report_message")
    def test_deletes_identical_and_warns_differing(self, mock_report_message):
        first_id, _identical_id = self._genotype_ids(self.identical_variant)
        differing_ids = self._genotype_ids(self.differing_variant)

        call_command("one_off_dedupe_cohort_genotypes", "--dry-run", stdout=StringIO())
        self.assertEqual(2, len(self._genotype_ids(self.identical_variant)))
        mock_report_message.assert_not_called()

        out = StringIO()
        call_command("one_off_dedupe_cohort_genotypes", stdout=out)
        self.assertEqual([first_id], self._genotype_ids(self.identical_variant))
        self.assertEqual(differing_ids, self._genotype_ids(self.differing_variant))
        self.assertIn(f"variant {self.differing_variant.pk}: rows {differing_ids}", out.getvalue())
        mock_report_message.assert_called_once()

        # "ERR" kept from the duplicated VCF, plus "EEO" from the other
        self.assertEqual((2, 3, 1), self._counts(self.identical_variant))
        self.assertEqual((9, 9, 9), self._counts(self.differing_variant), "No rows deleted, so not recounted")
