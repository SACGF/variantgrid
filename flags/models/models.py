import datetime
from collections import defaultdict
from functools import total_ordering, reduce
from operator import __and__
from typing import Tuple, List, Optional, Union, Dict, Iterable, Any, TypeVar

import django.dispatch
from django.contrib.auth.models import User
from django.db import models, transaction
from django.db.models import Q
from django.db.models.deletion import CASCADE, PROTECT
from django.db.models.expressions import Subquery
from django.db.models.query import QuerySet
from django.utils.timezone import now
from django_extensions.db.models import TimeStampedModel
from lazy import lazy

from flags.models.enums import FlagStatus
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.guardian_utils import admin_bot
from library.log_utils import report_message
from library.utils import empty_dict, ModelUtilsMixin, ChoicesEnum

flag_collection_extra_info_signal = django.dispatch.Signal()  # args: "flag_infos", "user"


@total_ordering
class FlagPermissionLevel(str, ChoicesEnum):
    NO_PERM = '0'
    USERS = 'U'
    OWNER = 'O'
    ADMIN = 'A'

    @property
    def level(self):
        if self == FlagPermissionLevel.USERS:
            return 1
        if self == FlagPermissionLevel.OWNER:
            return 2
        if self == FlagPermissionLevel.ADMIN:
            return 3
        return 0

    def __lt__(self, other):
        return self.level < other.level


class FlagResolution(TimeStampedModel, ModelUtilsMixin):
    id = models.TextField(primary_key=True)
    label = models.TextField()
    description = models.TextField()
    status = models.TextField(max_length=1, default=FlagStatus.OPEN, choices=FlagStatus.CHOICES)

    def __str__(self):
        return self.label


class FlagTypeContext(models.Model, ModelUtilsMixin):
    id = models.TextField(primary_key=True)
    label = models.TextField()

    def __str__(self):
        return self.label


class FlagType(TimeStampedModel, ModelUtilsMixin):
    id = models.TextField(primary_key=True)
    context = models.ForeignKey(FlagTypeContext, on_delete=CASCADE)
    label = models.TextField()
    description = models.TextField()
    # details on how to fix the flag (only shown to people who can edit it).
    help_text = models.TextField(default='')

    raise_permission = models.TextField(max_length=1, choices=FlagPermissionLevel.choices(), default=FlagPermissionLevel.ADMIN.value)

    def __str__(self):
        return self.label

    @property
    def raise_permission_enum(self) -> FlagPermissionLevel:
        return FlagPermissionLevel(self.raise_permission)

    permission = models.TextField(max_length=1, choices=FlagPermissionLevel.choices())

    @property
    def permission_enum(self) -> FlagPermissionLevel:
        return FlagPermissionLevel(self.permission)

    only_one = models.BooleanField(default=False)

    # can put days open until warning back in when we actually use it
    # days_open_until_warning = models.IntegerField(null=True, blank=True)
    attributes = models.JSONField(null=False, blank=True, default=empty_dict)
    comments_enabled = models.BooleanField(default=True)
    # determines on what screens the flag will be shown, more about if it causes clutter than importance
    importance = models.IntegerField(default=0)

    def resolution_for_status(self, status):
        ftr = FlagTypeResolution.objects.filter(
            flag_type=self,
            resolution__status=status
        ).first()
        if ftr:
            return ftr.resolution
        return None

    def default_resolution(self):
        res = self.resolution_for_status(FlagStatus.OPEN)
        if not res:
            res = self.resolution_for_status(FlagStatus.CLOSED)
        return res


class FlagTypeResolution(TimeStampedModel):
    flag_type = models.ForeignKey(FlagType, on_delete=CASCADE)
    resolution = models.ForeignKey(FlagResolution, on_delete=CASCADE)

    class Meta:
        unique_together = ('flag_type', 'resolution')


flag_comment_action = django.dispatch.Signal()  # args: "flag_comment", "old_resolution"


class Flag(TimeStampedModel):
    collection = models.ForeignKey('FlagCollection', on_delete=CASCADE)
    flag_type = models.ForeignKey(FlagType, on_delete=CASCADE)
    user_private = models.BooleanField(default=False)
    user = models.ForeignKey(User, on_delete=CASCADE, null=False)
    resolution = models.ForeignKey(FlagResolution, on_delete=PROTECT)
    data = models.JSONField(null=True, blank=True)

    @transaction.atomic()
    def flag_action(self,
                    resolution: FlagResolution = None,
                    user: User = None,
                    comment: str = None,
                    permission_check: bool = True,
                    first_comment: bool = False):

        if not user:
            permission_check = False
            user = admin_bot()

        current_permission = None
        if permission_check:
            current_permission = self.collection.permission_level(user)
            if self.user_private and (current_permission == FlagPermissionLevel.USERS and not self.user == user):
                raise PermissionError("User does not have permission to edit private flag")

        old_resolution = self.resolution
        if resolution == self.resolution and not first_comment:
            # only save resolutions that have changed
            resolution = None

        if resolution is not None:
            if permission_check:
                if current_permission < self.flag_type.permission_enum:
                    raise PermissionError(f"User does not have permission to {resolution.label} flag")

            # don't need to save the flag on creation as it will be created with the correct resolution
            if not first_comment:
                self.resolution = resolution
                self.save()

        if resolution is None and comment is None and not first_comment:
            # no-op
            return

        fc = FlagComment.objects.create(
            flag=self,
            user=user,
            text=comment or '',
            resolution=resolution
        )
        flag_comment_action.send(sender=Flag, flag_comment=fc, old_resolution=old_resolution)

    def __str__(self):
        return f"{self.flag_type} ({self.resolution})"


class FlagInfos:

    def __init__(self, flag_collections: List['FlagCollection'], flags: Iterable['Flag']):
        self.flag_collections = flag_collections
        self.flags_for_collection: Dict[Any, List[Flag]] = defaultdict(list)
        for flag in flags:
            self.flags_for_collection[flag.collection_id].append(flag)
        self._flagc_dict = {}
        self._flag_extra_info = {}
        for fc in self.flag_collections:
            fc._extra_info = {}
            self._flagc_dict[fc.pk] = fc
        self.sub_flag_types = []

    @lazy
    def ids(self):
        return [fc.pk for fc in self.flag_collections]

    def set_extra_info(self, pk, extra_info, source_object):
        fc = self._flagc_dict[pk]
        fc.extra_info = extra_info
        fc._source_object = source_object

    """
    # if we want to go back to providing extra info on the flag level
    def set_extra_flag_info(self, flag: Flag, extra_info: Dict[str,Any]):
        self._flag_extra_info[flag.id] = extra_info
    """

    def record_sub_flag(self, flag: Flag, label: str, letter: str):
        extra_info = self._flag_extra_info.get(flag.id, {'sub_flags':[]})
        self._flag_extra_info[flag.id] = extra_info
        sub_flags = extra_info['sub_flags']
        sub_flags.append({'label': label, 'letter': letter})

    def extra_flag_info(self, flag: Flag) -> Dict[str, Any]:
        return self._flag_extra_info.get(flag.id, {})


class FlagCollection(models.Model, GuardianPermissionsMixin):
    ADMIN_PERMISSION = 'admin_flagcollection'
    Q_OPEN_FLAGS = Q(resolution__status=FlagStatus.OPEN)
    REOPEN_IF_CLOSED_BY_BOT = 'reopen-if-bot'

    context = models.ForeignKey(FlagTypeContext, on_delete=PROTECT)

    class Meta:
        permissions = (
            # read/write permissions still exist by default
            ('admin_flagcollection', 'FlagCollection Admin'),
        )

    def __init__(self, *args, **kwargs):
        super(FlagCollection, self).__init__(*args, **kwargs)
        self._extra_info = None
        self._source_object = None

    @property
    def extra_info(self):
        if self._extra_info is None:
            fetch_flag_infos(flag_collections=[self], flags=[])
        return self._extra_info

    @extra_info.setter
    def extra_info(self, data: dict):
        self._extra_info = {**(self._extra_info or {}), **data}

    @property
    def source_object(self) -> Any:
        """
        The object that the FlagCollection is attached to, will be responsible for determining the user's permissions
        in relation to the FlagCollection
        """

        # ._source_object could be set either via getting FlagInfo (via a hook)
        # or by us directly going through the
        if not self._source_object:
            foreign_sets = [m for m in dir(self) if m.endswith('_set') and not m.startswith('flag')]
            for foreign_set in foreign_sets:
                source_object = getattr(self, foreign_set).first()
                if source_object:
                    self._source_object = source_object
                    break

        if not self._source_object:
            # appears to be an orphaned set
            report_message('Could not find source object for FlagCollection', extra_data={'flag_collection_id': self.id})

        return self._source_object

    def __str__(self):
        label = self.extra_info.get('label')
        if not label:
            label = f'FlagCollection object ({self.id})'
        return label

    def permission_level(self, user: User) -> FlagPermissionLevel:
        if user.is_superuser:
            # return owner instead of admin as admin currently really means automatic set by code
            return FlagPermissionLevel.OWNER

        so = self.source_object
        if not so:
            return FlagPermissionLevel.NO_PERM
        return so.flag_user_permission(user)

    def is_owner_or_admin(self, user: User) -> bool:
        permission = self.permission_level(user)
        return permission in (FlagPermissionLevel.ADMIN, FlagPermissionLevel.OWNER)

    def flags(self, user: User = None, only_open=False) -> QuerySet[Flag]:
        if not user:
            qs = Flag.objects.filter(collection=self)
        else:
            permission_level = self.permission_level(user)
            if not permission_level:
                return Flag.objects.none()
            qs = Flag.objects.filter(collection=self)
            if permission_level == FlagPermissionLevel.USERS:
                qs = qs.filter(Q(user_private=False) | Q(user=user))
        if only_open:
            qs = qs.filter(FlagCollection.Q_OPEN_FLAGS)
        return qs

    QST = TypeVar("QST", bound='FlagsMixin')

    @staticmethod
    def filter_for_open_flags(qs: QuerySet[QST], flag_types: Optional[List[FlagType]] = None) -> QuerySet[QST]:
        """
        @deprecated use filter_for_flags
        """
        return FlagCollection.filter_for_flags(qs=qs, flag_types=flag_types, open_only=True)

    @staticmethod
    def filter_for_flags(qs: QuerySet[QST], flag_types: Optional[List[FlagType]] = None, open_only: bool = True) -> QuerySet[QST]:
        """
        Takes the QuerySet and returns a filtered version where the item contain at least one of the provided flag_types
        e.g. if you passed in a QuerySet of Alleles and a missing 38 flag type, the resulting QuerySet will still be Alleles
        @param qs A QuerySet of models that have FlaxMixin
        @param flag_types if provided one of these flag types have to be included (otherwise any open flag can be included)
        @param open_only if True (default) closed flags will be ignored
        """
        flag_collections = Flag.objects.all()
        if open_only:
            flag_collections = flag_collections.filter(FlagCollection.Q_OPEN_FLAGS)
        if flag_types:
            flag_collections = flag_collections.filter(flag_type__in=flag_types)

        return qs.filter(flag_collection__in=Subquery(flag_collections.values('collection')))

    @staticmethod
    def filter_for_starred(qs: QuerySet, user: User) -> QuerySet:
        starred = FlagWatch.objects.filter(user=user).values('flag_collection').distinct()
        return qs.filter(flag_collection__in=starred)

    def close_open_flags_of_type(
            self,
            flag_type: FlagType,
            comment: str = None,
            user: User = None,
            resolution: FlagResolution = None,
            data: dict = None) -> int:
        """
        For admin usage, closes all open flags of a certain type
        @param flag_type close all flags of this type
        @param comment (optional) If provided, a FlagComment will be made for each flag we close
        @return: The number of flags closed, typically should be 0 or 1 but could be more
        """

        if not user:
            user = admin_bot()
        close_count = 0

        qs = Flag.objects.filter(collection=self, flag_type=flag_type).filter(FlagCollection.Q_OPEN_FLAGS)
        if data:
            for key, value in data.items():
                qs = qs.filter(**{f'data__{key}': value})

        for flag in qs:
            resolution = resolution if resolution else flag.flag_type.resolution_for_status(FlagStatus.CLOSED)
            flag.flag_action(user=user, comment=comment, resolution=resolution, permission_check=False)
            close_count = close_count + 1
        return close_count

    def get_open_flag_of_type(self, flag_type: FlagType = None, flag_type_attributes: dict = None):
        return self.get_flag_of_type(flag_type=flag_type, flag_type_attributes=flag_type_attributes, open_only=True)

    def get_flag_of_type(self, flag_type: FlagType = None, flag_type_attributes: dict = None, open_only=True) -> Optional[Flag]:
        if flag_type:
            qs = Flag.objects.filter(collection=self, flag_type=flag_type).order_by('-created')
            if open_only:
                qs = qs.filter(FlagCollection.Q_OPEN_FLAGS)
            return qs.first()
        if flag_type_attributes:
            ft_qs = FlagType.objects
            for key, value in flag_type_attributes.items():
                ft_qs = ft_qs.filter(**{f'attributes__{key}': value})
            qs = Flag.objects.filter(collection=self, flag_type__in=ft_qs).order_by('-created')
            if open_only:
                qs = qs.filter(FlagCollection.Q_OPEN_FLAGS)
            return qs.first()
        raise ValueError('Must provide flag_type or flag_type_attributes')

    def ensure_resolution(self, flag_type: Union[FlagType, str], resolution: Union[FlagResolution, str], comment: str = None):
        flag_type = FlagType.get(flag_type)
        resolution = FlagResolution.get(resolution)
        if not flag_type.only_one:
            raise ValueError('Can only call ensure_resolution on only_one flags')

        existing = Flag.objects.filter(collection=self, flag_type=flag_type).order_by('-created').first()
        if existing:
            if existing.resolution != resolution:
                existing.flag_action(resolution=resolution, comment=comment)
        else:
            if resolution.status == FlagStatus.OPEN:
                # don't add the flag if the resolution is Closed or Rejected
                self.add_flag(
                    resolution=resolution,
                    flag_type=flag_type,
                    comment=comment
                )

    @transaction.atomic()
    def get_or_create_open_flag_of_type(
            self,
            flag_type: FlagType,
            user: Optional[User] = None,
            comment: Optional[str] = None,
            user_private: bool = False,
            permission_check: bool = True,
            reopen: bool = False,
            reopen_if_bot_closed: bool = False,
            add_comment_if_open: bool = False,
            data: Optional[dict] = None,
            close_other_data: bool = False,
            only_if_new: bool = False) -> Tuple[Flag, bool]:
        """
        Returns the existing open flag or returns a new one
        :param flag_type: The type of flag to create
        :param user: The user creating the flag
        :param comment: A comment to be added onto the new flag
        :param user_private: (Not currently used)
        :param permission_check: Should we check to see if the user has permission to open the flag
        :param reopen: If True will re-open a closed flag (if there is one) rather than create a new flag. Will be treated as True for only_one types of flags.
        :param reopen_if_bot_closed: If True, and there's an existing flag that was closed by admin_bot, act as if reopen=True
        :param add_comment_if_open: If re-opening a closed flag, should the comment still be added?
        :param data: data for the flag, when looking for existing flags we check to see if they have this data.
        :param close_other_data: If we find a flag of the same type that has data different to the data provided, close it
        :param only_if_new: Only create a new flag if there isn't an existing one (in any state) for this type, data etc. Exclusive with reopen
        :return: A tuple of the flag and a boolean indicating if the flag was newly created
        """

        if flag_type.only_one and not only_if_new:
            reopen = True

        relevant_qs = Flag.objects.filter(collection=self, flag_type=flag_type).order_by('-created')

        if data:
            data_filters = []
            for key, value in data.items():
                data_filters.append(Q(**{f'data__{key}': value}))
            data_fitlers_q = reduce(__and__, data_filters)

            if close_other_data:
                close_us = relevant_qs.filter(FlagCollection.Q_OPEN_FLAGS).exclude(data_fitlers_q)
                for close_me in close_us:
                    resolution = close_me.flag_type.resolution_for_status(FlagStatus.CLOSED)
                    close_me.flag_action(user=user, comment='Data has changed. Raising a new flag.', resolution=resolution, permission_check=False)

            relevant_qs = relevant_qs.filter(data_fitlers_q)

        existing = relevant_qs.filter(FlagCollection.Q_OPEN_FLAGS).first()
        if existing:
            if add_comment_if_open:
                existing.flag_action(user=user, comment=comment)
            return existing, False

        existing = relevant_qs.first()
        if existing:
            if not reopen and reopen_if_bot_closed:
                # see if existing is closed, and was last closed by a bot
                if FlagComment.last(existing).user == admin_bot():
                    reopen = True

            if reopen:
                resolution = existing.flag_type.resolution_for_status(FlagStatus.OPEN)
                existing.flag_action(user=user, resolution=resolution, comment=comment, permission_check=permission_check)
                return existing, False
            if only_if_new:
                return existing, False

        return (self.add_flag(
            flag_type=flag_type,
            user=user,
            comment=comment,
            user_private=user_private,
            permission_check=permission_check,
            data=data), True)

    def add_flag(
            self,
            flag_type: FlagType,
            user: User = None,
            comment: str = None,
            user_private: bool = False,
            permission_check: bool = False,
            resolution: FlagResolution = None,
            data: dict = None) -> Flag:

        if not user:
            user = admin_bot()
            permission_check = False

        if permission_check:
            if flag_type.context_id != self.context_id:
                raise PermissionError(f"Flag type {flag_type.id} not available in flag context {self.context.label}")
            current_level = self.permission_level(user)
            required_level = flag_type.raise_permission
            if current_level < required_level:
                raise PermissionError(f"User does not have {required_level} permissions on flag collection")

        if resolution is None:
            resolution = flag_type.default_resolution()

        flag = Flag.objects.create(
            collection=self,
            flag_type=flag_type,
            user_private=user_private,
            resolution=resolution,
            user=user,
            data=data
        )
        flag.flag_action(comment=comment, user=user, resolution=resolution, permission_check=False, first_comment=True)

        return flag

    def set_watcher(self, user: User, watch: bool):
        if watch:
            FlagWatch.objects.get_or_create(flag_collection=self, user=user)
        else:
            FlagWatch.objects.filter(flag_collection=self, user=user).delete()

    def is_watcher(self, user: User):
        return FlagWatch.objects.filter(flag_collection=self, user=user).exists()

    def unseen_flag_activity(self, user: User) -> Optional[int]:
        """
        Returns None if not watching
        Otherwise returns how many flag comments that have not been seen (including 0)
        """
        fw = FlagWatch.objects.filter(flag_collection=self, user=user).first()
        if not fw:
            return None
        #FIXME filter to only flags that the user can see
        #excludes flags
        return FlagComment.objects.filter(flag__collection=self, created__gte=fw.since).exclude(user=user).count()


def fetch_flag_infos(flag_collections: List[FlagCollection], flags: Iterable[Flag], user: User = None) -> FlagInfos:
    flag_infos = FlagInfos(flag_collections=flag_collections, flags=flags)
    flag_collection_extra_info_signal.send(sender=FlagCollection, flag_infos=flag_infos, user=user)
    return flag_infos


class FlagWatch(models.Model):
    flag_collection = models.ForeignKey(FlagCollection, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)
    since = models.DateTimeField(default=now)

    class Meta:
        unique_together = ('flag_collection', 'user')


class FlagComment(TimeStampedModel):
    flag = models.ForeignKey(Flag, on_delete=CASCADE)
    user = models.ForeignKey(User, on_delete=CASCADE)
    text = models.TextField(blank=True)
    resolution = models.ForeignKey(FlagResolution, on_delete=PROTECT, null=True)

    @staticmethod
    def last(flag: Flag) -> 'FlagComment':
        return FlagComment.objects.filter(flag=flag).order_by('-created').first()


class FlagsMixin(models.Model):
    flag_collection = models.ForeignKey(FlagCollection, null=True, on_delete=CASCADE)

    class Meta:
        abstract = True

    def flag_type_context(self) -> FlagTypeContext:
        raise NotImplementedError('flag_type_context not implemented, required for FlagsMixin')

    def flag_user_permission(self, user: User) -> FlagPermissionLevel:
        if hasattr(self, 'can_write') and self.can_write(user):
            return FlagPermissionLevel.OWNER
        if hasattr(self, 'can_view'):
            return FlagPermissionLevel.USERS if self.can_view(user) else FlagPermissionLevel.NO_PERM
        return FlagPermissionLevel.USERS

    def flags_of_type(self, flag_type: FlagType) -> QuerySet[Flag]:
        if not self.flag_collection:
            return Flag.objects.none()
        return Flag.objects.filter(collection=self.flag_collection, flag_type=flag_type)

    @property
    def has_open_flags(self) -> bool:
        if not self.flag_collection:
            return False
        return Flag.objects.filter(collection=self.flag_collection).filter(FlagCollection.Q_OPEN_FLAGS).exists()

    def has_open_flag_with_attribute(self, attribute, value) -> bool:
        if not self.flag_collection:
            return False
        return Flag.objects.filter(collection=self.flag_collection).filter(FlagCollection.Q_OPEN_FLAGS) \
            .filter(**{'flag_type__attributes__%s' % attribute: value}).exists()

    def close_open_flags_of_type(
            self,
            flag_type: FlagType,
            comment=None,
            resolution: FlagResolution = None,
            data: dict = None):
        if not self.flag_collection:
            return False
        fc = self.flag_collection_safe
        fc.close_open_flags_of_type(flag_type, comment=comment, user=None, resolution=resolution, data=data)

    def has_flag_activity(self, since=None) -> bool:
        if since is None:
            since = datetime.datetime.now() - datetime.timedelta(days=1)
        if self.flag_collection:
            return FlagComment.objects.filter(flag__collection=self.flag_collection, created__gte=since).exists()
        return False

    @property
    def flag_collection_safe(self) -> FlagCollection:
        flag_collection = self.flag_collection
        if not flag_collection:
            flag_collection = FlagCollection.objects.create(context=self.flag_type_context())
            self.flag_collection = flag_collection
            self.save(update_fields=['flag_collection'])
        return flag_collection
