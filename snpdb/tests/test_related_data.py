from django.contrib.auth.models import User
from django.test import TestCase

from patients.models import Patient
from pedigree.models import (
    CohortSamplePedFileRecord,
    PedFile,
    PedFileFamily,
    PedFileRecord,
    Pedigree,
)
from snpdb.fake_data import (
    create_fake_cohort,
    create_fake_duo,
    create_fake_pedigree,
    create_fake_quad,
)
from snpdb.models import Cohort, Country, GenomeBuild, Lab, Organization, Trio
from snpdb.templatetags.related_data_tags import (
    related_data_for_cohort,
    related_data_for_patient,
    related_data_for_samples,
)


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


class RelatedDataForSamplesTest(TestCase):
    """ A sample's page lists the duos and quads it is a member of, as well as trios (#1889) """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='related_data_samples_user')[0]
        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.duo = create_fake_duo(cls.user, genome_build)
        cls.quad = create_fake_quad(cls.user, genome_build)

    def test_shows_duo_and_quad_with_member_sample(self):
        sibling = self.quad.sibling.sample
        related_data = related_data_for_samples({"user": self.user}, [self.duo.proband.sample, sibling])
        self.assertEqual(related_data["duos_and_samples"], [(self.duo, "proband")])
        self.assertEqual(related_data["quads_and_samples"], [(self.quad, "sibling")])


class RelatedDataPermissionsTest(TestCase):
    """ A patient's samples can come from another lab's VCF, and another lab can build its own sub cohort,
        trio and pedigree over this lab's samples - the tables only name what the user can view """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        organization = Organization.objects.get_or_create(name="Fake Org", group_name="fake_org")[0]
        country = Country.objects.get_or_create(name="Australia")[0]
        cohorts = []
        for lab_name in ["lab_a", "lab_b"]:
            lab = Lab.objects.get_or_create(name=lab_name, city="Adelaide", country=country,
                                            organization=organization, group_name=f"fake_org/{lab_name}")[0]
            user = User.objects.get_or_create(username=f"related_data_{lab_name}_user")[0]
            lab.group.user_set.add(user)
            cohorts.append(create_fake_cohort(user, cls.grch37, name=lab_name))
        cls.user, other_lab_user = [cohort.user for cohort in cohorts]
        cls.cohort, cls.other_lab_cohort = cohorts

        cls.patient = Patient.objects.create(first_name="RELATED", last_name="DATA")
        cls.sample = cls.cohort.cohortsample_set.get(sample__name="proband").sample
        cls.other_lab_sample = cls.other_lab_cohort.cohortsample_set.get(sample__name="proband").sample
        for sample in [cls.sample, cls.other_lab_sample]:
            sample.patient = cls.patient
            sample.save()

        cls.sub_cohort = cls.cohort.create_sub_cohort(cls.user, list(cls.cohort.get_samples()))
        cls.other_lab_sub_cohort = cls.cohort.create_sub_cohort(other_lab_user, list(cls.cohort.get_samples()))
        cls.other_lab_trio = cls._create_trio(other_lab_user, cls.other_lab_sub_cohort)
        cls.other_lab_pedigree = create_fake_pedigree(other_lab_user, cls.grch37, cohort=cls.other_lab_sub_cohort)
        for cohort_sample in cls.other_lab_sub_cohort.cohortsample_set.all():
            record = PedFileRecord.objects.create(family=cls.other_lab_pedigree.ped_file_family,
                                                  sample=cohort_sample.sample.name)
            CohortSamplePedFileRecord.objects.create(pedigree=cls.other_lab_pedigree, cohort_sample=cohort_sample,
                                                     ped_file_record=record)

    @staticmethod
    def _create_trio(user, cohort) -> Trio:
        cohort_samples = {cs.sample.name: cs for cs in cohort.cohortsample_set.all()}
        return Trio.objects.create(name=f"{cohort}_trio", user=user, cohort=cohort,
                                   mother=cohort_samples["mother"], mother_affected=True,
                                   father=cohort_samples["father"], father_affected=False,
                                   proband=cohort_samples["proband"])

    def test_patient_lists_only_visible_samples(self):
        related_data = related_data_for_patient({"user": self.user}, self.patient)
        self.assertEqual(list(related_data["samples"]), [self.sample])

    def test_samples_list_only_visible_cohorts(self):
        related_data = related_data_for_samples({"user": self.user}, self.patient.get_samples())
        self.assertEqual(list(related_data["samples"]), [self.sample])
        self.assertEqual([cohort for cohort, _ in related_data["cohorts_and_samples"]],
                         sorted([self.cohort, self.sub_cohort]))
        self.assertEqual(related_data["trios_and_samples"], [])
        self.assertEqual(related_data["pedigrees_and_samples"], [])

    def test_cohort_lists_only_visible_sub_cohorts(self):
        related_data = related_data_for_cohort({"user": self.user}, self.cohort)
        self.assertEqual(related_data["sub_cohorts"], [self.sub_cohort])
        self.assertEqual(related_data["trios"], [])
        self.assertEqual(related_data["pedigrees"], [])

    def test_cohort_filter_for_user_queryset(self):
        queryset = Cohort.objects.filter(pk__in=[self.cohort.pk, self.other_lab_cohort.pk, self.other_lab_sub_cohort.pk])
        self.assertEqual(list(Cohort.filter_for_user(self.user, queryset=queryset)), [self.cohort])
