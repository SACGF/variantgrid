from django.core.exceptions import ValidationError
from django.test import TestCase, override_settings

from classification.tests.models.test_utils import ClassificationTestUtils
from mme.matching import _our_patient_object
from mme.serializers.patient_profile import build_patient_profile
from mme.tests.fakes import (
    FakeClassification, FakeSubmission, FakeLab, FakeOrganization, make_term,
)
from ontology.models import OntologyService

SERVER_CONTACT = {"name": "Node Operator", "href": "mailto:node@example.org",
                  "institution": "Node Org"}


class LabMmeEnableValidationTestCase(TestCase):
    """ A lab cannot opt into MME without a resolvable contact (enable-time gate). """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, _ = ClassificationTestUtils.lab_and_user()

    def test_enable_without_contact_raises(self):
        self.lab.mme_enabled = True
        with self.assertRaises(ValidationError):
            self.lab.clean()

    def test_enable_with_contact_saves(self):
        self.lab.contact_email = "curator@lab.org"
        self.lab.mme_enabled = True
        self.lab.clean()   # no raise
        self.lab.save()
        self.lab.refresh_from_db()
        self.assertTrue(self.lab.mme_enabled)

    def test_disabled_without_contact_ok(self):
        self.lab.mme_enabled = False
        self.lab.clean()   # no raise


@override_settings(MME_ONTOLOGY_SNAKE_EXACT=True, MME_ONTOLOGY_PHENOTYPE_EXPANSION=False,
                   MME_ENABLED_PRODUCTION_SUBMIT=False)
class BuildProfileContactTestCase(TestCase):

    def _submission_with_hpo(self, lab):
        term = make_term("HP:0001250", OntologyService.HPO, 1250, "Seizure")
        classification = FakeClassification(terms=[term], lab=lab)
        return FakeSubmission(classification=classification)

    @override_settings(MME_CONTACT={"name": "", "href": ""})
    def test_uses_lab_contact(self):
        lab = FakeLab(name="Labby", contact_email="curator@lab.org",
                      organization=FakeOrganization(name="InstX"))
        profile = build_patient_profile(self._submission_with_hpo(lab))
        self.assertEqual(profile["contact"], {
            "name": "Labby",
            "href": "mailto:curator@lab.org",
            "institution": "InstX",
        })

    @override_settings(MME_CONTACT={"name": "", "href": ""})
    def test_raises_naming_lab_when_lab_and_server_empty(self):
        submission = self._submission_with_hpo(FakeLab(name="Labby"))   # no contact
        with self.assertRaises(ValueError) as ctx:
            build_patient_profile(submission)
        self.assertIn("Labby", str(ctx.exception))


@override_settings(MME_CONTACT=SERVER_CONTACT)
class OurPatientObjectContactTestCase(TestCase):
    """ Inbound serve: each classification is attributed to its own lab; a lab with no
        contact falls back to the server contact rather than dropping the match. """

    def test_distinct_lab_contacts(self):
        lab_a = FakeLab(name="Lab A", contact_email="a@lab.org",
                        organization=FakeOrganization(name="Org A"))
        lab_b = FakeLab(name="Lab B", contact_email="b@lab.org",
                        organization=FakeOrganization(name="Org B"))
        patient_a = _our_patient_object(FakeClassification(pk=1, lab=lab_a), None, [], [])
        patient_b = _our_patient_object(FakeClassification(pk=2, lab=lab_b), None, [], [])
        self.assertEqual(patient_a["contact"]["href"], "mailto:a@lab.org")
        self.assertEqual(patient_b["contact"]["href"], "mailto:b@lab.org")

    def test_falls_back_to_server_contact(self):
        classification = FakeClassification(lab=FakeLab(name="Labby"))   # no lab contact
        patient = _our_patient_object(classification, None, [], [])
        self.assertEqual(patient["contact"], SERVER_CONTACT)
