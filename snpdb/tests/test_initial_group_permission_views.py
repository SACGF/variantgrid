""" The initial group permission form on the organization page (views_lab.view_organization) """
from django.contrib.auth.models import User
from django.test import override_settings
from django.urls import reverse

from library.django_utils.unittest_utils import URLTestCase
from snpdb.models import (
    Lab,
    LabHead,
    Organization,
    OrganizationUserSettingsOverride,
    SettingsInitialGroupPermission,
)

ORG_GROUP_NAME = "test_initial_perm_org"


@override_settings(USER_SETTINGS_SHOW_GROUPS=True)
class OrganizationInitialGroupPermissionViewTest(URLTestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.organization = Organization.objects.create(name="Test Initial Perm Org", group_name=ORG_GROUP_NAME)
        cls.lab = Lab.objects.create(name="Test Initial Perm Lab", organization=cls.organization,
                                     group_name=f"{ORG_GROUP_NAME}/lab", city="Adelaide")
        cls.head = User.objects.create_user("initial_perm_org_head")
        cls.lab.group.user_set.add(cls.head)
        LabHead.objects.create(lab=cls.lab, user=cls.head)

    def test_org_head_saves_org_group_initial_permission(self):
        org_group = self.organization.group
        org_override = OrganizationUserSettingsOverride.objects.get_or_create(organization=self.organization)[0]

        self.client.force_login(self.head)
        url = reverse('view_organization', kwargs={"organization_id": self.organization.pk})
        prefix = f"{org_override.pk}_{org_group.name}"
        self.client.post(url, {f"{prefix}-read": "on", f"{prefix}-write": "on"})

        sigp = SettingsInitialGroupPermission.objects.get(settings=org_override, group=org_group)
        self.assertTrue(sigp.read)
        self.assertTrue(sigp.write)

