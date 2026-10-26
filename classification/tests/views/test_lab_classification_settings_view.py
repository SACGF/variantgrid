from django.contrib.auth.models import User
from django.urls import reverse

from library.django_utils.unittest_utils import URLTestCase, prevent_request_warnings
from snpdb.models import Lab, LabHead, Organization


class LabClassificationSettingsTabTest(URLTestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        organization = Organization.objects.create(name="Settings Tab Org", group_name="settings_tab_org")
        cls.lab = Lab.objects.create(name="Settings Tab Lab", organization=organization,
                                     group_name="settings_tab_org/settings_tab_lab", city="Adelaide")
        cls.head = User.objects.create_user("settings_tab_head")
        cls.member = User.objects.create_user("settings_tab_member")
        cls.outsider = User.objects.create_user("settings_tab_outsider")
        for user in (cls.head, cls.member):
            cls.lab.group.user_set.add(user)
        LabHead.objects.create(lab=cls.lab, user=cls.head)
        cls.url = reverse('lab_classification_settings_tab', kwargs={"pk": cls.lab.pk})

    def test_head_saves(self):
        self.client.force_login(self.head)
        response = self.client.post(self.url, {"consolidates_variant_classifications": "on"})
        self.assertEqual(response.status_code, 302)
        self.lab.classification_settings.refresh_from_db()
        self.assertTrue(self.lab.classification_settings.consolidates_variant_classifications)

    @prevent_request_warnings
    def test_access(self):
        for user, method, expected_status in [(self.member, "get", 200),
                                              (self.outsider, "get", 403),
                                              (self.member, "post", 403)]:
            with self.subTest(user=user.username, method=method):
                self.client.force_login(user)
                response = getattr(self.client, method)(self.url)
                self.assertEqual(response.status_code, expected_status)
