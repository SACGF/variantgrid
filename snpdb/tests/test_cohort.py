from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models import GenomeBuild, CohortGenotypeCollection, Cohort
from snpdb.tasks.cohort_genotype_tasks import cohort_genotype_task, create_cohort_genotype_collection, \
    create_cohort_genotype_and_launch_task
from snpdb.tests.test_data import create_fake_trio


class CohortGenotypeTestCase(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.trio1 = create_fake_trio(cls.user_owner, cls.grch37)
        cls.trio2 = create_fake_trio(cls.user_owner, cls.grch37)

    def test_sub_cohort(self):
        """ Create a cohort which should be a sub-cohort of trio1 """
        cohort = Cohort.objects.create(name="sub_cohort", genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        cohort.add_sample(samples[0].pk)
        cohort.add_sample(samples[1].pk)
        create_cohort_genotype_and_launch_task(cohort, run_async=False)
        self.assertTrue(cohort.is_sub_cohort())

    def test_not_sub_cohort(self):
        """ Create a cohort which is NOT a sub-cohort """
        cohort = Cohort.objects.create(name="not_a_sub_cohort", genome_build=self.grch37)
        samples = list(self.trio1.get_samples())
        cohort.add_sample(samples[0].pk)
        cohort.add_sample(samples[1].pk)
        create_cohort_genotype_and_launch_task(cohort, run_async=False)
        self.assertTrue(cohort.is_sub_cohort())  # Ok to be one

        samples2 = list(self.trio2.get_samples())
        cohort.add_sample(samples2[0].pk)  # Add a sample that is NOT in parent cohort, now can't be sub cohort
        self.assertFalse(cohort.is_sub_cohort())  # Shouldn't be a sub cohort

