from functools import lru_cache
from typing import Union, List

from django.conf import settings
from django.contrib.auth.models import Group, User
from django.core.exceptions import PermissionDenied
from django.db.models import Model, QuerySet
from guardian.shortcuts import get_groups_with_perms, get_users_with_perms, remove_perm, assign_perm


def is_superuser(user):
    """
    Can use this as the func in @user_passes_test annotation
    """
    return user.is_superuser


def all_users_group():
    g, _ = Group.objects.get_or_create(name=settings.LOGGED_IN_USERS_GROUP_NAME)
    return g


def public_group():
    g, _ = Group.objects.get_or_create(name=settings.PUBLIC_GROUP_NAME)
    return g


@lru_cache()
def admin_bot():
    return User.objects.get(username='admin_bot')


def bot_group():
    g, _ = Group.objects.get_or_create(name='variantgrid/bot')
    return g


def add_public_group_read_permission(obj):
    read_perm = DjangoPermission.perm(obj, DjangoPermission.READ)
    assign_perm(read_perm, public_group(), obj)


def _group_general_score(group_name):
    if group_name == settings.PUBLIC_GROUP_NAME:
        return 4
    if group_name == settings.LOGGED_IN_USERS_GROUP_NAME:
        return 3
    if '/' not in group_name:
        return 2
    return 1


def highest_group(group_names):
    highest_score = 0
    highest_group_name = None

    for name in group_names:
        score = _group_general_score(name)
        if score > highest_score:
            highest_score = score
            highest_group_name = name
    return highest_group_name


def groups_map(klass, ids, id_field='id'):
    id_to_group = {}
    if ids:
        if hasattr(ids[0], id_field):
            ids = [obj.id for obj in ids]

        objs = klass.objects.filter(pk__in=ids).only(id_field)
        for obj in objs:
            id_to_group[obj.id] = set(get_groups_with_perms(obj).values_list('name', flat=True))

    return id_to_group


def highest_group_map(klass, ids, id_field='id'):
    return [(key, highest_group(groups)) for key, groups in groups_map(klass, ids, id_field).items()]


class DjangoPermission:
    READ = 'view'
    WRITE = 'change'

    @staticmethod
    def perm(obj: Union[Model, List, QuerySet], permission):
        if isinstance(obj, QuerySet):
            klass = obj.model
        elif isinstance(obj, List):
            klass = obj[0]
        else:
            klass = obj

        return f"{permission}_{klass._meta.model_name}"


def assign_permission_to_user_and_groups(user: User, obj):
    """ adds to all non-public groups """

    from snpdb.models import UserSettings

    # Read permission to use and non public groups
    read_perm = DjangoPermission.perm(obj, DjangoPermission.READ)
    write_perm = DjangoPermission.perm(obj, DjangoPermission.WRITE)

    assign_perm(read_perm, user, obj)
    assign_perm(write_perm, user, obj)

    user_settings = UserSettings.get_for_user(user)
    read_groups, write_groups = user_settings.initial_perm_read_and_write_groups

    for group in read_groups:
        assign_perm(read_perm, group, obj)

    for group in write_groups:
        assign_perm(write_perm, group, obj)


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


def print_permissions(obj):
    """
    Prints all users and objects that have any permissions to the object
    """
    users = get_users_with_perms(obj)
    groups = get_groups_with_perms(obj)
    print(users)
    print(groups)


def check_can_write(obj, user):
    if not obj.can_write(user):
        raise PermissionDenied(f"You do not have WRITE permission for {obj.pk}")


def check_can_delete(user, pk, owner):
    can_delete = user and (user.is_superuser or user == owner)
    if not can_delete:
        raise PermissionDenied(f"You are not allowed to delete {pk}")
