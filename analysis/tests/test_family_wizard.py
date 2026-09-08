from types import SimpleNamespace

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from analysis.forms import UserDuoWizardForm, UserQuadWizardForm, UserTrioWizardForm
from analysis.models.enums import DuoSample, QuadSample, TrioSample
from analysis.views.views_wizard import _confident_sex
from library.django_utils.unittest_utils import URLTestCase
from patients.models_enums import Sex
from snpdb.models import Duo, DuoRelationship, GenomeBuild
from snpdb.tests.utils.fake_cohort_data import create_fake_duo


class ConfidentSexTest(TestCase):
    """ Which sex the wizard holds a sample's role to """

    @staticmethod
    def _sample(patient_sex: Sex, detected_sex: Sex):
        return SimpleNamespace(patient_sex=patient_sex, detected_sex=detected_sex)

    def test_patient_record_used_when_nothing_detected(self):
        self.assertEqual(Sex.FEMALE, _confident_sex(self._sample(Sex.FEMALE, Sex.UNKNOWN)))

    def test_falls_back_to_detected_sex(self):
        self.assertEqual(Sex.MALE, _confident_sex(self._sample(Sex.UNKNOWN, Sex.MALE)))

    def test_disagreement_leaves_every_role_open(self):
        self.assertEqual(Sex.UNKNOWN, _confident_sex(self._sample(Sex.FEMALE, Sex.MALE)))


class TrioWizardFormTest(TestCase):
    FEMALE_MALE_UNKNOWN = [Sex.FEMALE, Sex.MALE, Sex.UNKNOWN]

    def _form(self, data=None):
        return UserTrioWizardForm(data, sample_sexes=self.FEMALE_MALE_UNKNOWN)

    def test_sex_narrows_the_roles_on_offer(self):
        form = self._form()
        female_roles = [value for value, _ in form.fields["sample_1"].choices]
        male_roles = [value for value, _ in form.fields["sample_2"].choices]
        unknown_roles = [value for value, _ in form.fields["sample_3"].choices]

        self.assertEqual([TrioSample.MOTHER, TrioSample.PROBAND], female_roles)
        self.assertEqual([TrioSample.FATHER, TrioSample.PROBAND], male_roles)
        self.assertEqual([TrioSample.MOTHER, TrioSample.FATHER, TrioSample.PROBAND], unknown_roles)

    def test_female_cannot_be_posted_as_father(self):
        form = self._form({"sample_1": TrioSample.FATHER,
                           "sample_2": TrioSample.MOTHER,
                           "sample_3": TrioSample.PROBAND})
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
    def test_sibling_is_open_to_either_sex(self):
        form = UserQuadWizardForm(sample_sexes=[Sex.FEMALE, Sex.MALE, Sex.UNKNOWN, Sex.FEMALE])
        female_roles = [value for value, _ in form.fields["sample_1"].choices]
        self.assertEqual([QuadSample.MOTHER, QuadSample.PROBAND, QuadSample.SIBLING], female_roles)

    def test_affected_by_role_covers_the_sibling(self):
        form = UserQuadWizardForm({"sample_1": QuadSample.MOTHER,
                                   "sample_2": QuadSample.FATHER,
                                   "sample_3": QuadSample.PROBAND,
                                   "sample_4": QuadSample.SIBLING,
                                   "sample_4_affected": "on"},
                                  sample_sexes=[Sex.FEMALE, Sex.MALE, Sex.UNKNOWN, Sex.FEMALE])
        self.assertTrue(form.is_valid(), form.errors)
        self.assertTrue(form.affected_by_role[QuadSample.SIBLING])
        self.assertFalse(form.affected_by_role[QuadSample.MOTHER])


class DuoWizardFormTest(TestCase):
    """ A duo's roles are open to either sex - which parent it is comes from its own radio """

    def _form(self, data=None):
        return UserDuoWizardForm(data, sample_sexes=[Sex.FEMALE, Sex.MALE])

    def test_either_sample_can_take_either_role(self):
        form = self._form()
        for field in ("sample_1", "sample_2"):
            roles = [value for value, _ in form.fields[field].choices]
            self.assertEqual([DuoSample.PARENT, DuoSample.PROBAND], roles)

    def test_parent_relationship_and_affected_come_off_the_form(self):
        form = self._form({"sample_1": DuoSample.PARENT,
                           "sample_2": DuoSample.PROBAND,
                           "sample_1_affected": "on",
                           "relationship": DuoRelationship.FATHER})
        self.assertTrue(form.is_valid(), form.errors)
        self.assertEqual({DuoSample.PARENT: True, DuoSample.PROBAND: True}, form.affected_by_role)
        self.assertEqual(DuoRelationship.FATHER, form.parent_relationship)

    def test_two_samples_in_the_same_role_is_rejected(self):
        form = self._form({"sample_1": DuoSample.PROBAND,
                           "sample_2": DuoSample.PROBAND,
                           "relationship": DuoRelationship.MOTHER})
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
        cls.parent_sample_id = fixture_duo.parent.sample_id
        fixture_duo.delete()

    def test_post_creates_duo_from_roles(self):
        self.client.force_login(self.user)
        url = reverse("duo_wizard", kwargs={"cohort_id": self.cohort.pk,
                                            "sample1_id": self.proband_sample_id,
                                            "sample2_id": self.parent_sample_id})
        response = self.client.post(url, {"sample_1": DuoSample.PROBAND,
                                          "sample_2": DuoSample.PARENT,
                                          "sample_2_affected": "on",
                                          "relationship": DuoRelationship.FATHER,
                                          "proband_sex": ""})
        duo = Duo.objects.get(cohort=self.cohort)
        self.assertEqual(response.url, duo.get_absolute_url())
        self.assertEqual(duo.parent.sample_id, self.parent_sample_id)
        self.assertEqual(duo.proband.sample_id, self.proband_sample_id)
        self.assertEqual(duo.relationship, DuoRelationship.FATHER)
        self.assertTrue(duo.parent_affected)
