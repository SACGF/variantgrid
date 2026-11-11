from types import SimpleNamespace

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone

from analysis.enums import DuoSample, QuadSample, TrioSample
from analysis.forms import UserDuoWizardForm, UserQuadWizardForm, UserTrioWizardForm
from analysis.views.views_wizard import _confident_sex, family_name_suggestions
from library.django_utils.unittest_utils import URLTestCase
from patients.enums import Sex
from snpdb.enums import DuoRelationship
from snpdb.fake_data import create_fake_duo
from snpdb.models import VCF, Duo, GenomeBuild, Sample


class ConfidentSexTest(TestCase):
    """ Which sex the wizard holds a sample's role to """

    @staticmethod
    def _sample(patient_sex: Sex, detected_sex: Sex):
        return SimpleNamespace(patient_sex=patient_sex, detected_sex=detected_sex)

    def test_patient_record_used_when_nothing_detected(self):
        self.assertEqual(Sex.FEMALE, _confident_sex(self._sample(Sex.FEMALE, Sex.UNKNOWN)))

    def test_falls_back_to_detected_sex(self):
        self.assertEqual(Sex.MALE, _confident_sex(self._sample(Sex.UNKNOWN, Sex.MALE)))

    def test_patient_record_wins_a_disagreement(self):
        self.assertEqual(Sex.FEMALE, _confident_sex(self._sample(Sex.FEMALE, Sex.MALE)))


class FamilyNameSuggestionsTest(TestCase):
    """ #2760 - the VCF's name when it holds just the family, then what the sample names share """

    @classmethod
    def setUpTestData(cls):
        cls.user = User.objects.get_or_create(username='testuser_family_name')[0]

    def _samples(self, vcf_name: str, sample_names: list[str], extra_samples: int = 0) -> list[Sample]:
        vcf = VCF.objects.create(name=vcf_name, genotype_samples=len(sample_names) + extra_samples,
                                 genome_build=GenomeBuild.grch38(), user=self.user, date=timezone.now())
        samples = [Sample.objects.create(name=name, vcf=vcf) for name in sample_names]
        for i in range(extra_samples):
            Sample.objects.create(name=f"other_{i}", vcf=vcf)
        return samples

    def test_vcf_named_without_extension_then_common_sample_name(self):
        samples = self._samples("run7_trio.vcf.gz", ["FAM123_P", "FAM123_M", "FAM123_F"])
        self.assertEqual(["run7_trio", "FAM123"], family_name_suggestions(samples))

    def test_vcf_with_other_samples_is_not_suggested(self):
        samples = self._samples("run7_batch.vcf.gz", ["FAM123_P", "FAM123_M", "FAM123_F"], extra_samples=5)
        self.assertEqual(["FAM123"], family_name_suggestions(samples))

    def test_short_common_part_is_not_a_name(self):
        samples = self._samples("batch.vcf", ["AB_S1", "CD_S2", "EF_S3"], extra_samples=1)
        self.assertEqual([], family_name_suggestions(samples))


class TrioWizardFormTest(TestCase):
    def _form(self, data=None, patient_sexes=None):
        return UserTrioWizardForm(data, patient_sexes=patient_sexes)

    def test_patient_sex_withholds_the_other_parents_role(self):
        form = self._form(patient_sexes=[Sex.FEMALE, Sex.MALE, Sex.UNKNOWN])
        female_roles, male_roles, unknown_roles = ([value for value, _ in form.fields[f].choices]
                                                   for f in form.SAMPLE_FIELDS)
        self.assertEqual(["", TrioSample.MOTHER, TrioSample.PROBAND], female_roles)
        self.assertEqual(["", TrioSample.FATHER, TrioSample.PROBAND], male_roles)
        self.assertEqual(["", TrioSample.MOTHER, TrioSample.FATHER, TrioSample.PROBAND], unknown_roles)

    def test_female_patient_cannot_be_posted_as_father(self):
        form = self._form({"sample_1": TrioSample.FATHER,
                           "sample_2": TrioSample.MOTHER,
                           "sample_3": TrioSample.PROBAND},
                          patient_sexes=[Sex.FEMALE, Sex.UNKNOWN, Sex.UNKNOWN])
        self.assertFalse(form.is_valid())
        self.assertIn("sample_1", form.errors)

    def test_proband_is_affected_without_being_ticked(self):
        form = self._form({"sample_1": TrioSample.MOTHER,
                           "sample_2": TrioSample.FATHER,
                           "sample_3": TrioSample.PROBAND,
                           "sample_1_affected": "on"})
        self.assertTrue(form.is_valid(), form.errors)
        self.assertEqual({TrioSample.MOTHER: True,
                          TrioSample.FATHER: False,
                          TrioSample.PROBAND: True}, form.affected_by_role)
        self.assertEqual([TrioSample.MOTHER, TrioSample.FATHER, TrioSample.PROBAND], form.roles)

    def test_two_samples_in_the_same_role_is_rejected(self):
        form = self._form({"sample_1": TrioSample.MOTHER,
                           "sample_2": TrioSample.FATHER,
                           "sample_3": TrioSample.FATHER})
        self.assertFalse(form.is_valid())


class QuadWizardFormTest(TestCase):
    def test_affected_by_role_covers_the_sibling(self):
        form = UserQuadWizardForm({"sample_1": QuadSample.MOTHER,
                                   "sample_2": QuadSample.FATHER,
                                   "sample_3": QuadSample.PROBAND,
                                   "sample_4": QuadSample.SIBLING,
                                   "sample_4_affected": "on"})
        self.assertTrue(form.is_valid(), form.errors)
        self.assertTrue(form.affected_by_role[QuadSample.SIBLING])
        self.assertFalse(form.affected_by_role[QuadSample.MOTHER])


class DuoWizardFormTest(TestCase):
    """ The relative's role is which relative they are, so the Duo's relationship comes off the select """

    def _form(self, data=None):
        return UserDuoWizardForm(data)

    def test_relative_role_and_affected_come_off_the_form(self):
        form = self._form({"sample_1": DuoSample.FATHER,
                           "sample_2": DuoSample.PROBAND,
                           "sample_1_affected": "on"})
        self.assertTrue(form.is_valid(), form.errors)
        self.assertEqual({DuoSample.FATHER: True, DuoSample.PROBAND: True}, form.affected_by_role)
        self.assertEqual(DuoSample.FATHER, form.relative_role)

    def test_two_parents_without_a_proband_is_rejected(self):
        form = self._form({"sample_1": DuoSample.MOTHER,
                           "sample_2": DuoSample.FATHER})
        self.assertFalse(form.is_valid())

    def test_two_samples_in_the_same_role_is_rejected(self):
        form = self._form({"sample_1": DuoSample.PROBAND,
                           "sample_2": DuoSample.PROBAND})
        self.assertFalse(form.is_valid())


class DuoWizardViewTest(URLTestCase):
    """ The wizard hands back roles, not samples - check they land on the right side of the Duo """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username='testuser_duo_wizard')[0]
        # Borrow the fixture's cohort, then drop its Duo so the wizard makes its own
        fixture_duo = create_fake_duo(cls.user, GenomeBuild.get_name_or_alias("GRCh37"))
        cls.cohort = fixture_duo.cohort
        cls.proband_sample_id = fixture_duo.proband.sample_id
        cls.relative_sample_id = fixture_duo.relative.sample_id
        fixture_duo.delete()

    def _post_roles(self, relative_role, **extra):
        self.client.force_login(self.user)
        url = reverse("duo_wizard", kwargs={"cohort_id": self.cohort.pk,
                                            "sample1_id": self.proband_sample_id,
                                            "sample2_id": self.relative_sample_id})
        data = {"sample_1": DuoSample.PROBAND, "sample_2": relative_role, "proband_sex": "", **extra}
        return self.client.post(url, data)

    def test_post_creates_duo_from_roles(self):
        response = self._post_roles(DuoSample.FATHER, sample_2_affected="on")
        duo = Duo.objects.get(cohort=self.cohort)
        self.assertEqual(response.url, duo.get_absolute_url())
        self.assertEqual(duo.relative.sample_id, self.relative_sample_id)
        self.assertEqual(duo.proband.sample_id, self.proband_sample_id)
        self.assertEqual(duo.relationship, DuoRelationship.FATHER)
        self.assertTrue(duo.relative_affected)

    def test_post_creates_sibling_duo(self):
        self._post_roles(DuoSample.SIBLING, sample_2_affected="on")
        duo = Duo.objects.get(cohort=self.cohort)
        self.assertEqual(duo.relationship, DuoRelationship.SIBLING)
        self.assertTrue(duo.relative_is_sibling)
        self.assertTrue(duo.relative_affected)

    def test_post_uses_typed_name_or_falls_back_to_roles(self):
        self._post_roles(DuoSample.FATHER, name=" FAM123 ")
        self.assertEqual("FAM123", Duo.objects.get(cohort=self.cohort).name)
        Duo.objects.filter(cohort=self.cohort).delete()

        self._post_roles(DuoSample.FATHER, name="")
        # The fixture names its relative sample "mother", whatever role it is posted as
        self.assertEqual(f"mother/proband from {self.cohort}", Duo.objects.get(cohort=self.cohort).name)
