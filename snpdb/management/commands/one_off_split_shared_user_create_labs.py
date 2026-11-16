"""
Splits the shared lab of a settings.USER_CREATE_ORG_LABS holding organisation into a private lab per user.

Deployments used to put every self-registered user into one lab ("unknown/unknown") and the holding organisation's
group, so default sharing let each of them read the others' uploads. With a per-user lab pattern each user now gets
their own lab (snpdb/models/models_lab_request.py:create_user_org_labs). For each existing member this:
- adds them to their private lab and takes them out of the shared lab and the holding organisation's group
- re-grants each object permission the shared groups held to the private lab group of every member who has their own
  user permission on that object (uploads get both - library/guardian_utils.py:assign_permission_to_user_and_groups),
  then deletes the shared groups' permissions. Objects with no such user are reported: only an admin can see them now
- moves the shared lab's classifications to their creator's private lab, recomputing permissions and groupings
- points a default lab of the shared lab at the private lab

The emptied shared lab is left in place. --dry-run reports without writing.
Registered as a deploy-time ManualOperation (snpdb/migrations/0018_one_off_split_shared_user_create_labs.py).
"""
from collections import Counter

from django.conf import settings
from django.contrib.auth.models import Group, User
from django.contrib.contenttypes.models import ContentType
from django.core.management.base import BaseCommand
from django.db import transaction
from guardian.models import GroupObjectPermission, UserObjectPermission

from classification.models.classification import Classification, ClassificationModification
from classification.models.classification_grouping import ClassificationGrouping
from classification.services.overlaps_services import OverlapServices
from snpdb.models.models import Lab, Organization
from snpdb.models.models_lab_request import is_per_user_lab_pattern, user_create_org_lab
from snpdb.models.models_user_settings import UserSettingsOverride


class Command(BaseCommand):
    help = "Give each member of a shared USER_CREATE_ORG_LABS lab their own private lab, moving their data with them"
    category = "one-off"

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', help="Report what would change, change nothing")

    def handle(self, *args, **options):
        for org_group_name, lab_pattern in settings.USER_CREATE_ORG_LABS.items():
            if not is_per_user_lab_pattern(lab_pattern):
                continue
            if organization := Organization.objects.filter(group_name=org_group_name).first():
                with transaction.atomic():
                    self._split_organization(organization, lab_pattern, options["dry_run"])

    def _split_organization(self, organization: Organization, lab_pattern: str, dry_run: bool):
        private_labs: dict[User, Lab] = {}

        def private_lab_for(user: User) -> Lab:
            if user not in private_labs:
                private_labs[user] = user_create_org_lab(organization.group_name, lab_pattern, user)
            return private_labs[user]

        shared_labs = []
        for lab in organization.lab_set.all():
            members = list(lab.group.user_set.all())
            moving = [user for user in members if private_lab_for(user) != lab]
            if moving:
                shared_labs.append(lab)
        holding_group = organization.group
        holding_members = list(holding_group.user_set.all())
        for user in holding_members:
            private_lab_for(user)

        shared_groups = [lab.group for lab in shared_labs] + [holding_group]
        classification_qs = Classification.objects.filter(lab__in=shared_labs).select_related("user")
        for classification in classification_qs:
            private_lab_for(classification.user)
        shared_members = User.objects.filter(groups__in=[lab.group for lab in shared_labs]).distinct()

        self.stdout.write(f"{organization.group_name}: shared labs {', '.join(str(lab) for lab in shared_labs) or 'none'}, "
                          f"{shared_members.count()} shared lab members, {len(holding_members)} organisation group "
                          f"members, {classification_qs.count()} classifications to move")
        new_permissions, orphaned = self._regranted_permissions(shared_groups, private_labs)
        self.stdout.write(f"{len(new_permissions)} permissions to re-grant to private labs")
        for content_type, count in orphaned.most_common():
            self.stdout.write(f"  {count} {content_type} permissions with no owner to re-grant to (admin-only after this)")
        if dry_run:
            transaction.set_rollback(True)  # private_lab_for created labs
            return

        for user, private_lab in private_labs.items():
            private_lab.group.user_set.add(user)
        GroupObjectPermission.objects.bulk_create(new_permissions, ignore_conflicts=True, batch_size=5000)
        GroupObjectPermission.objects.filter(group__in=shared_groups).delete()

        with OverlapServices.discordance_notifications_suppressed(), OverlapServices.overlap_recalcs_deferred():
            for classification in classification_qs:
                private_lab = private_lab_for(classification.user)
                Classification.objects.filter(pk=classification.pk).update(lab=private_lab)
                classification.lab = private_lab
                classification.fix_permissions(fix_modifications=True)
                ClassificationGrouping.assign_grouping_for_classification(classification)
            for grouping in ClassificationGrouping.objects.filter(dirty=True, lab__organization=organization):
                grouping.update()

        for user, private_lab in private_labs.items():
            UserSettingsOverride.objects.filter(user=user, default_lab__in=shared_labs).update(default_lab=private_lab)
        for group in shared_groups:
            group.user_set.clear()
        self.stdout.write(f"Moved {len(private_labs)} users into private labs")

    @staticmethod
    def _regranted_permissions(shared_groups: list[Group],
                               private_labs: dict[User, Lab]) -> tuple[list[GroupObjectPermission], Counter]:
        """ The shared groups' object permissions for each private lab whose user has the same permission themselves,
            and per content type, how many have no such user. Classifications are left to fix_permissions """
        classification_types = ContentType.objects.get_for_models(Classification, ClassificationModification).values()
        shared_qs = GroupObjectPermission.objects.filter(group__in=shared_groups).exclude(
            content_type__in=classification_types)
        shared = set(shared_qs.values_list("content_type_id", "object_pk", "permission_id"))

        private_group_ids = {user.pk: lab.group.pk for user, lab in private_labs.items()}
        user_permissions_qs = UserObjectPermission.objects.filter(user_id__in=private_group_ids).values_list(
            "user_id", "content_type_id", "object_pk", "permission_id")
        new_permissions = []
        regranted = set()
        for user_id, content_type_id, object_pk, permission_id in user_permissions_qs.iterator():
            key = (content_type_id, object_pk, permission_id)
            if key in shared:
                regranted.add(key)
                new_permissions.append(GroupObjectPermission(group_id=private_group_ids[user_id],
                                                             content_type_id=content_type_id,
                                                             object_pk=object_pk, permission_id=permission_id))

        content_types = ContentType.objects.in_bulk()
        orphaned = Counter(str(content_types[content_type_id]) for content_type_id, _, _ in shared - regranted)
        return new_permissions, orphaned
