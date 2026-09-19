from typing import Optional

from django.contrib.auth.models import User
from django.db import connection
from django.test import TestCase
from django.test.utils import CaptureQueriesContext

from annotation.fake_annotation import get_fake_annotation_version
from library.django_utils.django_partition import temporary_db_table
from snpdb.models import (
    CohortGenotype,
    CohortGenotypeCollection,
    CohortGenotypeCollectionType,
    CohortGenotypeCommonFilterVersion,
    CommonVariantClassified,
    GenomeBuild,
)
from snpdb.tasks.cohort_genotype_tasks import (
    cohort_genotype_task,
    common_variant_classified_task,
    create_cohort_genotype_collection,
)
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort, create_fake_trio
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class CohortGenotypeTestCase(TestCase):

    def test_create_cohort(self):
        """ Recreates issue #307 - @see https://github.com/SACGF/variantgrid/issues/307 """

        user_owner = User.objects.get_or_create(username='testuser')[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        trio = create_fake_trio(user_owner, grch37)
        cohort = trio.cohort
        CohortGenotypeCollection.objects.filter(cohort=cohort, cohort_version=cohort.version).delete()
        cgc = create_cohort_genotype_collection(cohort)
        cohort_genotype_task(cgc.pk)

    def test_zygosity_filter_deleted_sample(self):
        """ Recreates issue #613 - @see https://github.com/SACGF/variantgrid/issues/613

            When generating a regex, we used to not take into account deleted samples, thus the regex
            could be less than the sample_zygosity length, and could 'shift' and match wrong things
        """

        user_owner = User.objects.get_or_create(username='testuser')[0]
        grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(grch37)
        trio = create_fake_trio(user_owner, grch37)
        cohort = trio.cohort
        CohortGenotypeCollection.objects.filter(cohort=cohort, cohort_version=cohort.version).delete()
        cgc = create_cohort_genotype_collection(cohort)
        cohort_genotype_task(cgc.pk)

        # Delete a sample
        cgc.cohort.vcf.sample_set.first().delete()
        # Ensure that the number of "." wildcard entries match number of ORIGINAL samples
        zygosity = cgc.get_sample_zygosity_regex({}, {})
        self.assertEqual(len(zygosity), cgc.num_samples)


class SampleGenotypeCopyNumberTest(TestCase):
    """ The caller's copy number / copy ratio, which has no packed column and is read back out of
        the CohortGenotype JSON """

    def setUp(self):
        user = User.objects.get_or_create(username='testuser')[0]
        self.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        self.cohort = create_fake_cohort(user, self.genome_build)
        self.sample = self.cohort.cohortsample_set.get(sample__name="proband").sample
        self.variant = slowly_create_test_variant("3", 128198980, 'G', 'C', self.genome_build)

    def _sample_genotype(self, copy_number_field: Optional[str], sample_format: list,
                         info: Optional[dict] = None):
        vcf = self.sample.vcf
        vcf.copy_number_field = copy_number_field
        vcf.save()
        cg = CohortGenotype.objects.create(collection=self.cohort.cohort_genotype_collection,
                                           variant=self.variant, samples_zygosity="OOO",
                                           format=sample_format, info=info or {})
        return cg.get_sample_genotype(self.sample)

    def test_reads_the_declared_field_from_this_sample_s_format(self):
        sample_genotype = self._sample_genotype("SM", [{"SM": [4.31428]}, {"SM": [1.0]}, {"SM": [1.1]}])
        self.assertEqual(4.31428, sample_genotype.copy_number_value)

    def test_falls_back_to_info_for_a_single_sample_vcf(self):
        sample_genotype = self._sample_genotype("CN", [{}, {}, {}], info={"CN": 12})
        self.assertEqual(12, sample_genotype.copy_number_value)

    def test_is_none_when_the_vcf_declares_no_field(self):
        sample_genotype = self._sample_genotype(None, [{"SM": [4.31428]}, {}, {}])
        self.assertIsNone(sample_genotype.copy_number_value)

    def test_is_none_when_the_record_carries_no_value(self):
        sample_genotype = self._sample_genotype("SM", [{"BC": [22]}, {}, {}])
        self.assertIsNone(sample_genotype.copy_number_value)


class CommonVariantClassifiedTest(TestCase):
    """ Classifying a common variant moves its CohortGenotype rows out of the common partition into the
        uncommon one, so they stop being skipped by the common filter optimisation. """

    def setUp(self):
        user = User.objects.get_or_create(username='testuser')[0]
        self.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(self.genome_build)
        self.cohort = create_fake_cohort(user, self.genome_build)
        self.uncommon_cgc = self.cohort.cohort_genotype_collection

        self.common_filter = CohortGenotypeCommonFilterVersion.objects.create(
            gnomad_version="2.1.1", gnomad_af_min=0.05, genome_build=self.genome_build)
        self.common_cgc = CohortGenotypeCollection.objects.create(
            cohort=self.cohort, cohort_version=self.cohort.version,
            num_samples=self.uncommon_cgc.num_samples,
            collection_type=CohortGenotypeCollectionType.COMMON,
            common_filter=self.common_filter)
        self.uncommon_cgc.common_collection = self.common_cgc
        self.uncommon_cgc.save()

        self.variant = slowly_create_test_variant("3", 128198980, 'G', 'C', self.genome_build)
        with temporary_db_table(CohortGenotype, self.common_cgc.get_partition_table()):
            CohortGenotype.objects.create(collection=self.common_cgc, variant=self.variant,
                                          samples_zygosity="OOO", het_count=0, hom_count=3)

    def _partition_count(self, cgc) -> int:
        with temporary_db_table(CohortGenotype, cgc.get_partition_table()):
            return CohortGenotype.objects.filter(collection=cgc).count()

    def test_moves_genotypes_from_the_common_to_the_uncommon_partition(self):
        self.assertEqual(self._partition_count(self.common_cgc), 1)
        self.assertEqual(self._partition_count(self.uncommon_cgc), 0)

        common_variant_classified_task(self.variant.pk, self.common_filter.pk)

        self.assertEqual(self._partition_count(self.common_cgc), 0)
        self.assertEqual(self._partition_count(self.uncommon_cgc), 1)
        self.assertTrue(CommonVariantClassified.objects.filter(variant=self.variant,
                                                               common_filter=self.common_filter).exists())

    def test_delete_does_not_expand_over_the_inheritance_parent(self):
        """ A delete off snpdb_cohortgenotype locks every other collection's partition, which deadlocks
            concurrent VCF imports (SACGF/variantgrid_sapath#450) """
        with CaptureQueriesContext(connection) as queries:
            common_variant_classified_task(self.variant.pk, self.common_filter.pk)

        base_table = '"snpdb_cohortgenotype"'
        tree_wide = [q["sql"] for q in queries
                     if q["sql"].lstrip().upper().startswith("DELETE") and base_table in q["sql"]]
        self.assertEqual(tree_wide, [], f"Delete expanded over the inheritance tree: {tree_wide}")
