from io import StringIO

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
        cls.variant = slowly_create_test_variant("3", 128198980, 'A', 'T', grch37)
        cls.collection = VariantZygosityCountCollection.get_global_germline_counts()
        create_variant_zygosity_counts()

        duplicated_cohort = create_fake_cohort(user, grch37, name="dedupe_duplicated")
        cls.cgc = duplicated_cohort.cohort_genotype_collection
        # Like the partitions that hold duplicates: created without the unique (collection, variant) index
        partition_table = cls.cgc.get_partition_table()
        run_sql(f'ALTER TABLE "{partition_table}" DROP CONSTRAINT "{partition_table}_collection_id_variant_id_key"')
        with temporary_db_table(CohortGenotype, partition_table):
            make_cohort_genotype(cls.cgc, cls.variant, "ERR")
            make_cohort_genotype(cls.cgc, cls.variant, "OOE")
        other_cohort = create_fake_cohort(user, grch37, name="dedupe_other")
        make_cohort_genotype(other_cohort.cohort_genotype_collection, cls.variant, "EEO")

        for cohort in [duplicated_cohort, other_cohort]:
            update_variant_zygosity_count_for_vcf(cls.collection, cohort.vcf, '+')
        # The count UPDATE joins whichever duplicate Postgres picks - here the later "OOE" one
        VariantZygosityCount.objects.filter(collection=cls.collection, variant=cls.variant).update(
            ref_count=0, het_count=3, hom_count=3)

    def _genotype_ids(self) -> list[int]:
        return list(CohortGenotype.objects.filter(collection=self.cgc).order_by("pk").values_list("pk", flat=True))

    def test_keeps_first_row_and_recounts(self):
        first_id, _later_id = self._genotype_ids()

        call_command("dedupe_cohort_genotypes", "--dry-run", stdout=StringIO())
        self.assertEqual(2, len(self._genotype_ids()))

        out = StringIO()
        call_command("dedupe_cohort_genotypes", stdout=out)
        self.assertEqual([first_id], self._genotype_ids())
        self.assertIn("1 with differing genotypes", out.getvalue())

        vzc = VariantZygosityCount.objects.get(collection=self.collection, variant=self.variant)
        # "ERR" kept from the duplicated VCF, plus "EEO" from the other
        self.assertEqual((2, 3, 1), (vzc.ref_count, vzc.het_count, vzc.hom_count))
