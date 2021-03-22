from django.contrib.auth.models import User
from django.test import TestCase

from snpdb.models import GenomeBuild, CohortGenotypeCollection
from snpdb.tasks.cohort_genotype_tasks import cohort_genotype_task, create_cohort_genotype_collection
from snpdb.tests.test_data import create_fake_trio


class AlleleTestCase(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.user_owner = User.objects.get_or_create(username='testuser')[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.trio = create_fake_trio(cls.user_owner, cls.grch37)

    def test_create_cohort(self):
        """ Recreates issue #307 - @see https://github.com/SACGF/variantgrid/issues/307 """
        cohort = self.trio.cohort
        CohortGenotypeCollection.objects.filter(cohort=cohort, cohort_version=cohort.version).delete()
        cgc = create_cohort_genotype_collection(cohort)
        cohort_genotype_task(cgc.pk)
