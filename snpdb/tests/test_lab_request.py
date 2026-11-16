from io import StringIO

from django.contrib.auth.models import Group, User
from django.core import signing
from django.core.management import call_command
from django.test import TestCase, override_settings
from django.urls import reverse

from classification.enums import ShareLevel, SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from classification.tests.models.test_utils import ClassificationTestUtils
from library.request_context import set_thread_variable
from snpdb.enums import LabRequestStatus
from snpdb.fake_data import create_fake_cohort
from snpdb.lab_request_approval import approve_lab_request
from snpdb.models import GenomeBuild, Lab, LabHead, Organization
from snpdb.models.models_lab_request import (
    LabRequest,
    create_user_org_labs,
    email_matches_website,
    get_user_private_lab,
)

HOLDING = "test_holding"
PRIVATE_LABS = {HOLDING: "user_group_%(username)s"}


def _create_user(username: str, **kwargs) -> User:
    """ The user post_save handler isn't connected in unit tests """
    user = User.objects.create_user(username, **kwargs)
    create_user_org_labs(user)
    return user


class EmailMatchesWebsiteTest(TestCase):

    def test_email_matches_website(self):
        for email, website, matches in [
            ("someone@rsu.lv", "https://www.rsu.lv", True),
            ("someone@RSU.lv", "https://rsu.lv/about", True),
            ("someone@mail.rsu.lv", "https://www.rsu.lv", True),
            ("someone@rsu.lv", "https://genetics.rsu.lv", True),
            ("someone@notrsu.lv", "https://rsu.lv", False),
            ("someone@gmail.com", "https://rsu.lv", False),
            ("someone@rsu.lv", "https://rsu.lv.example.com", False),
        ]:
            self.assertEqual(email_matches_website(email, website), matches, f"{email} / {website}")


@override_settings(USER_CREATE_ORG_LABS=PRIVATE_LABS, USER_CREATE_ORG_MESSAGE={})
class PrivateLabTest(TestCase):

    def setUp(self):
        set_thread_variable('request', None)
        self.alice = _create_user("lab_request_alice", email="alice@rsu.lv")
        self.bob = _create_user("lab_request_bob", email="bob@rsu.lv")

    def test_new_users_get_private_labs(self):
        alice_lab = get_user_private_lab(self.alice)
        self.assertNotEqual(alice_lab, get_user_private_lab(self.bob))
        self.assertEqual(list(Lab.valid_labs_qs(self.alice)), [alice_lab])
        self.assertFalse(self.alice.groups.filter(name=HOLDING).exists(), "Doesn't join the holding org's group")

        cohort = create_fake_cohort(self.alice, GenomeBuild.grch37())
        self.assertTrue(cohort.vcf.can_view(self.alice))
        self.assertFalse(cohort.vcf.can_view(self.bob))

    def test_approve_renames_private_lab(self):
        private_lab = get_user_private_lab(self.alice)
        group_pk = private_lab.group.pk
        cohort = create_fake_cohort(self.alice, GenomeBuild.grch37())
        ClassificationTestUtils.setUp()
        data = {SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}, SpecialEKeys.ALLELE_ORIGIN: {"value": "germline"}}
        classification = Classification.create(user=self.alice, lab=private_lab, data=data, save=True,
                                                source=SubmissionSource.API)
        classification.publish_latest(user=self.alice, share_level=ShareLevel.ORGANISATION)

        lab_request = LabRequest.objects.create(user=self.alice, organization_name="Riga Uni", lab_name="Genetics",
                                                website="https://www.rsu.lv", email="alice@rsu.lv")
        self.assertEqual(LabRequest.from_confirmation_token(lab_request.confirmation_token), lab_request)
        with self.assertRaises(signing.BadSignature):
            LabRequest.from_confirmation_token(lab_request.confirmation_token + "x")
        lab_request.confirm_email()

        admin = User.objects.create_superuser("lab_request_admin")
        with self.captureOnCommitCallbacks(execute=True):
            lab = approve_lab_request(lab_request, admin, None, "Riga Uni", "Genetics")

        self.assertEqual(lab.pk, private_lab.pk)
        self.assertEqual(lab.group_name, "riga_uni/genetics")
        self.assertEqual(lab.group.pk, group_pk, "Same group row, so its object permissions stand")
        self.assertIsNone(get_user_private_lab(self.alice))
        self.assertTrue(LabHead.objects.filter(lab=lab, user=self.alice).exists())
        self.assertTrue(self.alice.groups.filter(name="riga_uni").exists())
        self.assertTrue(cohort.vcf.can_view(lab.group))
        self.assertEqual(LabRequest.objects.get(pk=lab_request.pk).status, LabRequestStatus.APPROVED)

        organization_group = Organization.objects.get(group_name="riga_uni").group
        modification = classification.last_published_version
        self.assertTrue(modification.can_view(organization_group), "Organisation share follows the new org")
        self.assertFalse(modification.can_view(Group.objects.get(name=HOLDING)))

        self.assertIn(self.bob, lab.candidate_members_qs(self.alice), "Colleague on the request's email domain")


class SplitSharedUserCreateLabTest(TestCase):

    def test_split_shared_lab(self):
        set_thread_variable('request', None)
        with override_settings(USER_CREATE_ORG_LABS={HOLDING: "shared"}, USER_CREATE_ORG_MESSAGE={}):
            alice = _create_user("lab_request_alice")
            bob = _create_user("lab_request_bob")
            cohort = create_fake_cohort(alice, GenomeBuild.grch37())
        shared_lab = Lab.objects.get(group_name=f"{HOLDING}/shared")
        self.assertTrue(cohort.vcf.can_view(bob), "Shared lab - everyone reads everyone's uploads")

        with override_settings(USER_CREATE_ORG_LABS=PRIVATE_LABS):
            call_command("one_off_split_shared_user_create_labs", "--dry-run", stdout=StringIO())
            self.assertTrue(cohort.vcf.can_view(bob), "Dry run changes nothing")
            self.assertIsNone(get_user_private_lab(alice))

            call_command("one_off_split_shared_user_create_labs", stdout=StringIO())
            alice_lab = get_user_private_lab(alice)
            self.assertTrue(cohort.vcf.can_view(alice_lab.group))
            self.assertFalse(cohort.vcf.can_view(User.objects.get(pk=bob.pk)))
            self.assertFalse(shared_lab.group.user_set.exists())
            self.assertFalse(Group.objects.get(name=HOLDING).user_set.exists())
            self.assertEqual(list(Lab.valid_labs_qs(bob)), [get_user_private_lab(bob)])


@override_settings(USER_CREATE_ORG_LABS=PRIVATE_LABS, USER_CREATE_ORG_MESSAGE={})
class LabRequestViewsTest(TestCase):

    def test_request_confirm_approve(self):
        set_thread_variable('request', None)
        alice = _create_user("lab_request_alice", email="alice@rsu.lv")
        admin = User.objects.create_superuser("lab_request_admin")
        self.client.force_login(alice)
        self.assertContains(self.client.get(reverse("view_user_settings")), reverse("lab_request"))

        post = {"organization_name": "Riga Uni", "lab_name": "Genetics", "website": "https://www.rsu.lv",
                "email": "alice@gmail.com"}
        self.assertContains(self.client.post(reverse("lab_request"), post), "same domain")
        post["email"] = "alice@rsu.lv"
        self.client.post(reverse("lab_request"), post)
        lab_request = LabRequest.objects.get(user=alice)
        self.assertContains(self.client.get(reverse("lab_request")), "confirm your email")

        self.client.get(reverse("lab_request_confirm", kwargs={"token": lab_request.confirmation_token}))
        self.assertContains(self.client.get(reverse("lab_request")), "waiting for an administrator")

        self.client.force_login(admin)
        self.assertContains(self.client.get(reverse("lab_requests")), "Riga Uni")
        url = reverse("view_lab_request", kwargs={"pk": lab_request.pk})
        self.assertContains(self.client.get(url), "lab-request-approve-form")
        response = self.client.post(url, {"action": "approve", "organization_name": "Riga Uni", "lab_name": "Genetics"})
        lab = Lab.objects.get(group_name="riga_uni/genetics")
        self.assertRedirects(response, reverse("view_lab", kwargs={"lab_id": lab.pk}), fetch_redirect_response=False)
