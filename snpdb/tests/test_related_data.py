from django.contrib.auth.models import User
from django.test import TestCase

from pedigree.models import PedFile, PedFileFamily, Pedigree
from snpdb.models import GenomeBuild, Trio
from snpdb.templatetags.related_data_tags import related_data_for_cohort
from snpdb.tests.utils.fake_cohort_data import create_fake_cohort


class RelatedDataForCohortTest(TestCase):
    """ Sub cohorts are made off the VCF page, then trios/pedigrees made off those - the parent
        cohort page needs to show them too """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='related_data_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.cohort = create_fake_cohort(cls.user, genome_build)
        cls.sub_cohort = cls.cohort.create_sub_cohort(cls.user, list(cls.cohort.get_samples()))

        cohort_samples = {cs.sample.name: cs for cs in cls.sub_cohort.cohortsample_set.all()}
        cls.trio = Trio.objects.create(name="sub_cohort_trio", user=cls.user, cohort=cls.sub_cohort,
                                       mother=cohort_samples["mother"], mother_affected=True,
                                       father=cohort_samples["father"], father_affected=False,
                                       proband=cohort_samples["proband"])
        ped_file = PedFile.objects.create(name="sub_cohort_ped", user=cls.user)
        ped_file_family = PedFileFamily.objects.create(ped_file=ped_file, name="family")
        cls.pedigree = Pedigree.objects.create(name="sub_cohort_pedigree", user=cls.user,
                                               cohort=cls.sub_cohort, ped_file_family=ped_file_family)

    def _related_data(self):
        return related_data_for_cohort({"user": self.user}, self.cohort)

    def test_shows_sub_cohorts(self):
        self.assertEqual(self._related_data()["sub_cohorts"], [self.sub_cohort])

    def test_shows_sub_cohort_trios(self):
        self.assertEqual(self._related_data()["trios"], [self.trio])

    def test_shows_sub_cohort_pedigrees(self):
        self.assertEqual(self._related_data()["pedigrees"], [self.pedigree])
