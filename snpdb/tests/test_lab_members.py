from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.test import TestCase, override_settings
from django.urls import reverse
from threadlocals.threadlocals import set_thread_variable

from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from snpdb.models import Lab, LabHead, Organization, UserSettingsOverride
from user_messages.models import Message

ORG_GROUP_NAME = "test_lab_members_org"
LAB_GROUP_NAME = f"{ORG_GROUP_NAME}/test_lab_members_lab"
OTHER_LAB_GROUP_NAME = f"{ORG_GROUP_NAME}/test_lab_members_other_lab"
ORG_MESSAGE = {ORG_GROUP_NAME: "A lab head will add you to a lab"}


def _create_lab_data():
    """ An org with 2 labs - a head/member/outsider in the org, and someone with no org at all """
    organization = Organization.objects.create(name="Test Members Org", group_name=ORG_GROUP_NAME)
    lab = Lab.objects.create(name="Test Members Lab", organization=organization,
                             group_name=LAB_GROUP_NAME, city="Adelaide")
    other_lab = Lab.objects.create(name="Test Members Other Lab", organization=organization,
                                   group_name=OTHER_LAB_GROUP_NAME, city="Adelaide")
    return organization, lab, other_lab


class LabMembershipModelTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.organization, cls.lab, cls.other_lab = _create_lab_data()

        cls.head = User.objects.create_user("lab_members_head")
        cls.member = User.objects.create_user("lab_members_member")
        cls.org_user = User.objects.create_user("lab_members_org_user")
        cls.inactive_org_user = User.objects.create_user("lab_members_inactive", is_active=False)
        cls.outsider = User.objects.create_user("lab_members_outsider")
        cls.superuser = User.objects.create_superuser("lab_members_superuser")

        for user in (cls.head, cls.member, cls.org_user, cls.inactive_org_user):
            cls.organization.group.user_set.add(user)
        for user in (cls.head, cls.member):
            cls.lab.group.user_set.add(user)
        LabHead.objects.create(lab=cls.lab, user=cls.head)

    def setUp(self):
        # AdminNotificationBuilder logs an Event against the thread-local user - a client request in
        # an earlier test leaves a rolled-back one behind
        set_thread_variable('request', None)
        self.lab.refresh_from_db()

    @override_settings(USER_CREATE_ORG_MESSAGE=ORG_MESSAGE)
    def test_organization_add_member_grants_no_lab(self):
        user = User.objects.create_user("lab_members_new_org_user")
        self.organization.add_member(user)

        self.assertTrue(user.groups.filter(name=ORG_GROUP_NAME).exists())
        self.assertFalse(Lab.valid_labs_qs(user).exists())
        self.assertEqual(Message.objects.filter(recipient=user).count(), 1)

        self.organization.add_member(user)
        self.assertEqual(Message.objects.filter(recipient=user).count(), 1, "2nd add doesn't re-message")

    def test_can_manage_members_setting_off(self):
        self.assertFalse(self.lab.can_manage_members(self.head))
        self.assertFalse(self.lab.can_manage_members(self.superuser))

    @override_settings(LAB_HEAD_MANAGE_MEMBERS=True)
    def test_can_manage_members_setting_on(self):
        self.assertTrue(self.lab.can_manage_members(self.head))
        self.assertTrue(self.lab.can_manage_members(self.superuser))
        self.assertFalse(self.lab.can_manage_members(self.member))
        self.assertFalse(self.other_lab.can_manage_members(self.head))

    def test_candidate_members_qs(self):
        candidates = set(self.lab.candidate_members_qs(self.head))
        self.assertEqual(candidates, {self.org_user}, "Org users not already in the lab")

        superuser_candidates = set(self.lab.candidate_members_qs(self.superuser))
        self.assertIn(self.outsider, superuser_candidates, "Superuser can add from outside the org")
        self.assertNotIn(self.member, superuser_candidates, "Existing members excluded")
        self.assertNotIn(self.inactive_org_user, superuser_candidates, "Inactive users excluded")

    def test_add_member(self):
        self.lab.add_member(self.org_user, added_by=self.head)

        self.assertTrue(self.org_user.groups.filter(name=LAB_GROUP_NAME).exists())
        user_settings_override = UserSettingsOverride.objects.get(user=self.org_user)
        self.assertEqual(user_settings_override.default_lab, self.lab)

    def test_remove_member_keeps_org_and_drops_head(self):
        LabHead.objects.create(lab=self.lab, user=self.member)
        self.lab.remove_member(self.member, removed_by=self.head)

        self.assertFalse(self.member.groups.filter(name=LAB_GROUP_NAME).exists())
        self.assertTrue(self.member.groups.filter(name=ORG_GROUP_NAME).exists(), "Org membership stands")
        self.assertFalse(LabHead.objects.filter(lab=self.lab, user=self.member).exists())

    def test_head_cannot_remove_last_head_or_self(self):
        with self.assertRaises(PermissionDenied):
            self.lab.remove_member(self.head, removed_by=self.member)
        with self.assertRaises(PermissionDenied):
            self.lab.remove_member(self.head, removed_by=self.head)

        self.lab.remove_member(self.head, removed_by=self.superuser)
        self.assertFalse(self.head.groups.filter(name=LAB_GROUP_NAME).exists())

    def test_external_membership(self):
        with override_settings(LAB_EXTERNAL_MEMBERSHIP={LAB_GROUP_NAME: "Ask IT"}):
            self.assertEqual(self.lab.external_membership, "Ask IT")
            with self.assertRaises(PermissionDenied):
                self.lab.remove_member(self.member, removed_by=self.head)

            self.lab.add_member(self.org_user, added_by=self.head)
            self.assertTrue(self.org_user.groups.filter(name=LAB_GROUP_NAME).exists(), "Adding still works")

            self.lab.remove_member(self.member, removed_by=self.superuser)
            self.assertFalse(self.member.groups.filter(name=LAB_GROUP_NAME).exists())


@override_settings(LAB_HEAD_MANAGE_MEMBERS=True)
class LabMembershipViewTest(URLTestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.organization, cls.lab, cls.other_lab = _create_lab_data()

        cls.head = User.objects.create_user("lab_members_view_head")
        cls.member = User.objects.create_user("lab_members_view_member")
        cls.org_user = User.objects.create_user("lab_members_view_org_user")
        cls.other_lab_head = User.objects.create_user("lab_members_view_other_head")
        cls.superuser = User.objects.create_superuser("lab_members_view_superuser")

        for user in (cls.head, cls.member, cls.org_user, cls.other_lab_head):
            cls.organization.group.user_set.add(user)
        for user in (cls.head, cls.member):
            cls.lab.group.user_set.add(user)
        LabHead.objects.create(lab=cls.lab, user=cls.head)
        LabHead.objects.create(lab=cls.other_lab, user=cls.other_lab_head)

    def test_head_can_add_and_remove(self):
        self.client.force_login(self.head)

        response = self.client.post(reverse('lab_add_member', kwargs={"pk": self.lab.pk}),
                                    {"user": self.org_user.pk})
        self.assertEqual(response.status_code, 302)
        self.assertTrue(self.org_user.groups.filter(name=LAB_GROUP_NAME).exists())

        response = self.client.post(reverse('lab_remove_member', kwargs={"pk": self.lab.pk}),
                                    {"user": self.org_user.pk})
        self.assertEqual(response.status_code, 302)
        self.assertFalse(self.org_user.groups.filter(name=LAB_GROUP_NAME).exists())

    def test_view_lab_renders(self):
        self.client.force_login(self.head)
        response = self.client.get(reverse('view_lab', kwargs={"lab_id": self.lab.pk}))
        self.assertEqual(response.status_code, 200)
        self.assertTrue(response.context["can_manage_members"])

    def test_members_tab(self):
        # lab_members_tab is in URLS_NAME_REGISTER, which is fixed when settings load, so only a
        # superuser gets past that wrapper here - the head's access is covered by can_manage_members
        self.client.force_login(self.superuser)
        response = self.client.get(reverse('lab_members_tab', kwargs={"pk": self.lab.pk}))
        self.assertEqual(response.status_code, 200)

    @prevent_request_warnings
    def test_non_head_denied(self):
        for user in (self.member, self.other_lab_head):
            self.client.force_login(user)
            with self.subTest(user=user.username):
                response = self.client.post(reverse('lab_add_member', kwargs={"pk": self.lab.pk}),
                                            {"user": self.org_user.pk})
                self.assertEqual(response.status_code, 403)

    @prevent_request_warnings
    def test_autocomplete_requires_head(self):
        url = reverse('lab_add_member_autocomplete') + f"?forward=%7B%22lab_id%22%3A%20{self.lab.pk}%7D"
        self.client.force_login(self.member)
        self.assertEqual(self.client.get(url).status_code, 403)

        self.client.force_login(self.head)
        self.assertEqual(self.client.get(url).status_code, 200)


class LabMembersUrlsDisabledTest(URLTestCase):
    """ With LAB_HEAD_MANAGE_MEMBERS off the feature doesn't exist for anyone """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.organization, cls.lab, _ = _create_lab_data()
        cls.head = User.objects.create_user("lab_members_disabled_head")
        cls.lab.group.user_set.add(cls.head)
        LabHead.objects.create(lab=cls.lab, user=cls.head)

    @prevent_request_warnings
    def test_views_denied(self):
        self.client.force_login(self.head)
        for url in (reverse('lab_members_tab', kwargs={"pk": self.lab.pk}),
                    reverse('lab_add_member', kwargs={"pk": self.lab.pk}),
                    reverse('lab_remove_member', kwargs={"pk": self.lab.pk})):
            with self.subTest(url=url):
                self.assertEqual(self.client.post(url).status_code, 403)
