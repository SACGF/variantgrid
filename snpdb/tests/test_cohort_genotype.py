from django.contrib.auth.models import User
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.models import GenomeBuild, CohortGenotypeCollection
from snpdb.tasks.cohort_genotype_tasks import cohort_genotype_task, create_cohort_genotype_collection
from snpdb.tests.test_data import create_fake_trio


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
