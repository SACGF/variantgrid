from django.test import TestCase, override_settings

from mme.contact import (
    lab_mme_contact,
    settings_mme_contact,
    mme_contact_for_classification,
)
from mme.tests.fakes import FakeClassification, FakeLab, FakeOrganization

SERVER_CONTACT = {"name": "Node Operator", "href": "mailto:node@example.org",
                  "institution": "Node Org"}


class LabMmeContactTestCase(TestCase):

    def test_full_lab_contact(self):
        lab = FakeLab(name="Labby", contact_name="Dr Curator",
                      contact_email="curator@lab.org", url="https://lab.org",
                      organization=FakeOrganization(name="InstX"))
        self.assertEqual(lab_mme_contact(lab), {
            "name": "Dr Curator",
            "href": "mailto:curator@lab.org",
            "institution": "InstX",
        })

    def test_falls_back_to_lab_name_and_url(self):
        # No contact_name -> lab.name; no email but a URL -> href from URL.
        lab = FakeLab(name="Labby", url="https://lab.org",
                      organization=FakeOrganization(name="InstX"))
        self.assertEqual(lab_mme_contact(lab), {
            "name": "Labby",
            "href": "https://lab.org",
            "institution": "InstX",
        })

    def test_institution_falls_back_to_lab_name(self):
        lab = FakeLab(name="Labby", contact_email="curator@lab.org")
        contact = lab_mme_contact(lab)
        self.assertEqual(contact["institution"], "Labby")

    def test_no_href_returns_empty(self):
        # name resolvable (lab.name) but neither email nor a URL -> {}.
        lab = FakeLab(name="Labby", url="not-a-url")
        self.assertEqual(lab_mme_contact(lab), {})

    def test_nothing_set_returns_empty(self):
        self.assertEqual(lab_mme_contact(FakeLab()), {})

    def test_none_lab_returns_empty(self):
        self.assertEqual(lab_mme_contact(None), {})


class SettingsMmeContactTestCase(TestCase):

    @override_settings(MME_CONTACT=SERVER_CONTACT)
    def test_valid_settings_contact(self):
        self.assertEqual(settings_mme_contact(), SERVER_CONTACT)

    @override_settings(MME_CONTACT={"name": "", "href": ""})
    def test_incomplete_settings_contact_empty(self):
        self.assertEqual(settings_mme_contact(), {})


class MmeContactForClassificationTestCase(TestCase):

    @override_settings(MME_CONTACT=SERVER_CONTACT)
    def test_lab_contact_wins(self):
        lab = FakeLab(name="Labby", contact_email="curator@lab.org",
                      organization=FakeOrganization(name="InstX"))
        classification = FakeClassification(lab=lab)
        self.assertEqual(mme_contact_for_classification(classification), {
            "name": "Labby",
            "href": "mailto:curator@lab.org",
            "institution": "InstX",
        })

    @override_settings(MME_CONTACT=SERVER_CONTACT)
    def test_falls_back_to_server_when_lab_empty(self):
        classification = FakeClassification(lab=FakeLab(name="Labby"))
        self.assertEqual(mme_contact_for_classification(classification), SERVER_CONTACT)

    @override_settings(MME_CONTACT={"name": "", "href": ""})
    def test_empty_when_both_empty(self):
        classification = FakeClassification(lab=FakeLab(name="Labby"))
        self.assertEqual(mme_contact_for_classification(classification), {})
