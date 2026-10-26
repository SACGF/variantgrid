"""
django-guardian plumbing: DjangoPermission builds the `view_<model>` / `change_<model>` permission
strings, the standard groups (all_users_group, public_group, bot_group) and the `admin_bot` system
user are looked up here, user_in_group / user_in_group_or_superuser are the group-membership checks, and assign_permission_to_user_and_groups grants a new object's read/write
perms from the user's UserSettings initial groups (the one sanctioned import of snpdb from library).
Models check permissions through library/django_utils/guardian_permissions_mixin.py, not directly here.
"""
from functools import lru_cache
from typing import Union

from django.conf import settings
from django.contrib.auth.models import Group, Permission, User
from django.core.exceptions import PermissionDenied
from django.db.models import Model, QuerySet
from guardian.ctypes import get_content_type
from guardian.models import GroupObjectPermission, UserObjectPermission
from guardian.shortcuts import assign_perm, get_groups_with_perms, get_users_with_perms, remove_perm


def is_superuser(user):
    """
    Can use this as the func in @user_passes_test annotation
    """
    return user.is_superuser


def user_in_group(user: User, group: Union[Group, str]) -> bool:
    group_name = group.name if isinstance(group, Group) else group
    return user.groups.filter(name=group_name).exists()


def user_in_group_or_superuser(user: User, group: Union[Group, str]) -> bool:
    return user.is_superuser or user_in_group(user, group)


def all_users_group():
    g, _ = Group.objects.get_or_create(name=settings.LOGGED_IN_USERS_GROUP_NAME)
    return g


def public_group():
    g, _ = Group.objects.get_or_create(name=settings.PUBLIC_GROUP_NAME)
    return g


@lru_cache
def _cached_admin_bot() -> User:
    return User.objects.get(username='admin_bot')


def admin_bot() -> User:
    """
    lru_cache is great for performance, but it can break unit tests by caching a User object across sessions (when the DB has been rolled back).
    Specifically this object is often referenced by things that want to save the user to the DB.
    """
    if settings.UNIT_TEST:
        return User.objects.get(username='admin_bot')
    else:
        return _cached_admin_bot()


def bot_group() -> Group:
    g, _ = Group.objects.get_or_create(name='variantgrid/bot')
    return g


def add_public_group_read_permission(obj):
    read_perm = DjangoPermission.perm(obj, DjangoPermission.READ)
    assign_perm(read_perm, public_group(), obj)


class DjangoPermission:
    READ = 'view'
    WRITE = 'change'

    @staticmethod
    def perm(obj: Union[Model, list, QuerySet], permission):
        if isinstance(obj, QuerySet):
            klass = obj.model
        elif isinstance(obj, list):
            klass = obj[0]
        else:
            klass = obj

        return f"{permission}_{klass._meta.model_name}"


def assign_permission_to_user_and_groups(user: User, obj):
    """ adds to all non-public groups

        One Permission query and a bulk_create per permission: guardian's assign_perm resolves the
        ContentType and Permission and does a get_or_create per row, which a VCF import pays for
        every sample. ignore_conflicts leans on the same unique constraint get_or_create did. """

    from snpdb.models import UserSettings

    # Read permission to user and non-public groups
    ctype = get_content_type(obj)
    codenames = [DjangoPermission.perm(obj, DjangoPermission.READ),
                 DjangoPermission.perm(obj, DjangoPermission.WRITE)]
    permissions = {p.codename: p for p in Permission.objects.filter(content_type=ctype, codename__in=codenames)}
    missing = set(codenames) - permissions.keys()
    if missing:
        raise Permission.DoesNotExist(f"No {ctype.app_label} permission(s) {sorted(missing)}")
    read_perm, write_perm = (permissions[codename] for codename in codenames)

    user_settings = UserSettings.get_for_user(user)
    read_groups, write_groups = user_settings.initial_perm_read_and_write_groups

    for permission, groups in [(read_perm, read_groups), (write_perm, write_groups)]:
        UserObjectPermission.objects.assign_perm_to_many(permission, [user], obj, ignore_conflicts=True)
        if groups:
            GroupObjectPermission.objects.assign_perm_to_many(permission, groups, obj, ignore_conflicts=True)


def clear_permissions(obj, permissions):
    """
    Removes the given permissions from the object for any user and any object
    that might have them
    """
    users = get_users_with_perms(obj, with_group_users=False)
    for user in users:
        for permission in permissions:
            remove_perm(permission, user, obj)

    groups = get_groups_with_perms(obj)
    for group in groups:
        for permission in permissions:
            remove_perm(permission, group, obj)


def check_can_write(obj, user: User):
    if not obj.can_write(user):
        raise PermissionDenied(f"You do not have WRITE permission for {obj.pk}")


def check_can_delete(user: User, pk, owner):
    can_delete = user and (user.is_superuser or user == owner)
    if not can_delete:
        raise PermissionDenied(f"You are not allowed to delete {pk}")
