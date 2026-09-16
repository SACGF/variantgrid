from typing import Optional

from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import CohortGenotype, CohortGenotypeCollection, GenomeBuild
from snpdb.tasks.cohort_genotype_tasks import (
    cohort_genotype_task,
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
