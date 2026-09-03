from types import SimpleNamespace

from django.test import TestCase

from analysis.forms import UserQuadWizardForm, UserTrioWizardForm
from analysis.models.enums import QuadSample, TrioSample
from analysis.views.views_wizard import _confident_sex
from patients.models_enums import Sex


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
