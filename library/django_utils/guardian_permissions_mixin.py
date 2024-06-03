import logging
from typing import Union

from django.conf import settings
from django.contrib.auth.models import Group, User
from django.core.exceptions import PermissionDenied
from django.shortcuts import get_object_or_404
from guardian.shortcuts import get_objects_for_user, get_objects_for_group, get_group_perms

from library.guardian_utils import assign_permission_to_user_and_groups, DjangoPermission


class GuardianPermissionsMixin:

    @classmethod
    def get_read_perm(cls):
        perm_cls = cls.get_permission_class()
        return DjangoPermission.perm(perm_cls, DjangoPermission.READ)

    @classmethod
    def get_write_perm(cls):
        perm_cls = cls.get_permission_class()
        return DjangoPermission.perm(perm_cls, DjangoPermission.WRITE)

    @classmethod
    def get_permission_class(cls):
        """ Object can use another classes permissions """
        return cls

    def get_permission_object(self):
        """ Object can use another objects permissions """
        return self

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        """ Object can use another objects permissions """
        return queryset

    def can_view(self, user_or_group: Union[User, Group]) -> bool:
        if not user_or_group:
            return False
        perm_obj = self.get_permission_object()
        if perm_obj != self:
            return perm_obj.can_view(user_or_group)

        if isinstance(user_or_group, Group):
            return self.get_read_perm() in get_group_perms(user_or_group, self)
        return user_or_group.has_perm(self.get_read_perm(), self)

    def can_write(self, user_or_group: Union[User, Group]) -> bool:
        if not user_or_group:
            return False
        perm_obj = self.get_permission_object()
        if perm_obj != self:
            return perm_obj.can_write(user_or_group)
        if isinstance(user_or_group, Group):
            return self.get_write_perm() in get_group_perms(user_or_group, self)
        return user_or_group.has_perm(self.get_write_perm(), self)

    def check_can_write(self, user_or_group: Union[User, Group]):
        if not self.can_write(user_or_group):
            msg = f"You do not have WRITE permission for {self.pk}"
            raise PermissionDenied(msg)

    def check_can_view(self, user_or_group: Union[User, Group]):
        if not self.can_view(user_or_group):
            msg = f"You do not have READ permission to view {self.pk}"
            raise PermissionDenied(msg)

    @classmethod
    def filter_for_user(cls, user, queryset=None, **kwargs):
        # QuerySet evaluates to False if it has no values, so check against None specifically
        klass = cls.get_permission_class()
        if klass != cls:
            # logging.info("%s delegating to %s", cls, klass)
            queryset = klass.filter_for_user(user, queryset=queryset, **kwargs)
        else:
            if queryset is not None:
                klass = queryset

            if user and user.is_authenticated:
                queryset = get_objects_for_user(user, cls.get_read_perm(), klass=klass, accept_global_perms=True)
            else:
                # No user - try public (non-logged in users) access
                group = Group.objects.get(name=settings.PUBLIC_GROUP_NAME)
                queryset = get_objects_for_group(group, cls.get_read_perm(), klass=klass, accept_global_perms=True)

        return cls._filter_from_permission_object_qs(queryset)

    @classmethod
    def get_for_user(cls, user, pk, write=False):
        obj = get_object_or_404(cls, pk=pk)
        if write:
            obj.check_can_write(user)
        else:
            obj.check_can_view(user)
        return obj


class GuardianPermissionsAutoInitialSaveMixin(GuardianPermissionsMixin):
    """ Automatically assigns permissions on initial save (unless you specify 'assign_permissions')
        This *must* be inherited from before Model in the class definition to call this save not model's """

    def save(self, **kwargs):
        assign_permissions = kwargs.pop("assign_permissions", None)
        initial_save = not self.pk
        super().save(**kwargs)
        if assign_permissions is None:
            assign_permissions = initial_save
        if assign_permissions:
            if user := getattr(self, "user", None):
                assign_permission_to_user_and_groups(user, self)
            else:
                raise ValueError(f"{self} tried to set permissions without a user")
