"""
assign_permission_to_user_and_groups writes the object permission rows every new user-owned object
starts with: the owner's read and write, plus the read/write groups from their UserSettings initial
group permissions. The rows are written in bulk, so what lands is worth pinning.
"""
from django.contrib.auth.models import Group, User
from django.contrib.contenttypes.models import ContentType
from django.test import TestCase
from guardian.models import GroupObjectPermission, UserObjectPermission

from analysis.models import Analysis
from library.guardian_utils import assign_permission_to_user_and_groups
from snpdb.models import GenomeBuild, GlobalSettings, SettingsInitialGroupPermission


class AssignPermissionToUserAndGroupsTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="guardian_utils_user")
        # A Group post_save handler gives each new group GlobalSettings read=True, write=False
        cls.read_group = Group.objects.create(name="guardian_utils_read")
        cls.read_write_group = Group.objects.create(name="guardian_utils_read_write")
        cls.user.groups.add(cls.read_group, cls.read_write_group)
        SettingsInitialGroupPermission.objects.filter(settings=GlobalSettings.objects.get(),
                                                      group=cls.read_write_group).update(write=True)
        # Analysis is a GuardianPermissionsAutoInitialSaveMixin, so creating it assigns
        cls.analysis = Analysis.objects.create(name="guardian utils", user=cls.user,
                                               genome_build=GenomeBuild.grch37())

    def _rows(self, model, holder_field) -> set[tuple[str, str]]:
        content_type = ContentType.objects.get_for_model(Analysis)
        qs = model.objects.filter(content_type=content_type, object_pk=str(self.analysis.pk))
        return {(str(getattr(row, holder_field)), row.permission.codename) for row in qs}

    def test_user_and_initial_groups_get_their_rows(self):
        self.assertEqual(self._rows(UserObjectPermission, "user"),
                         {("guardian_utils_user", "view_analysis"),
                          ("guardian_utils_user", "change_analysis")})
        self.assertEqual(self._rows(GroupObjectPermission, "group"),
                         {("guardian_utils_read", "view_analysis"),
                          ("guardian_utils_read_write", "view_analysis"),
                          ("guardian_utils_read_write", "change_analysis")})

    def test_assigning_again_does_not_duplicate(self):
        """ bulk_create(ignore_conflicts) in place of get_or_create - a re-assign has to stay a no-op """
        user_rows = self._rows(UserObjectPermission, "user")
        group_rows = self._rows(GroupObjectPermission, "group")
        assign_permission_to_user_and_groups(self.user, self.analysis)
        self.assertEqual(self._rows(UserObjectPermission, "user"), user_rows)
        self.assertEqual(self._rows(GroupObjectPermission, "group"), group_rows)
        content_type = ContentType.objects.get_for_model(Analysis)
        self.assertEqual(UserObjectPermission.objects.filter(content_type=content_type,
                                                             object_pk=str(self.analysis.pk)).count(), 2)
