import datetime
from typing import Iterable, Dict, Any, Union, List, Optional

from django.contrib.auth.models import User
from django.utils import timezone
from django.utils.timezone import now
from rest_framework.response import Response
from rest_framework.views import APIView

from flags.models import Flag, FlagComment, FlagType, FlagCollection, FlagWatch
from flags.models.enums import FlagStatus
from flags.models.models import FlagResolution, FlagTypeResolution, fetch_flag_infos
from library.django_utils import ensure_timezone_aware
from library.utils import empty_to_none
from snpdb.models import Lab
from snpdb.models.models_user_settings import UserSettings


class CommentDetails:

    def __init__(self, helper, comment: FlagComment):
        self.comment = comment
        self.helper = helper
        self.detailed = False

    def to_json(self) -> Dict[str, Any]:
        status = None
        resolution = self.comment.resolution
        if resolution:
            status = resolution.status
            resolution = resolution.id

        json_entry = {
            'id': self.comment.id,
            'user': self.comment.user_id,
            'flag': self.comment.flag_id,
            'created': self.comment.created.timestamp(),
            'status': status,  # FIXME, don't transmit status on comments
            'resolution': resolution
        }
        if self.detailed:
            json_entry['text'] = self.comment.text
        else:
            summary = self.comment.text
            if len(summary) < 60:
                json_entry['text'] = summary
            else:
                json_entry['summary'] = summary[0:60] + '...'
        return json_entry


class FlagHelper:

    def __init__(self, flag_collections: Union[FlagCollection, Iterable[FlagCollection]], user: User):

        if flag_collections is None:
            flag_collections = []
        elif not isinstance(flag_collections, Iterable):
            flag_collections = [flag_collections]
        elif not isinstance(flag_collections, list):
            flag_collections = list(flag_collections)

        self.since = None
        self.user = user
        self.flag_collections: List[FlagCollection] = flag_collections
        self.flag_types: Dict[Any, FlagType] = {}
        self.flag_resolutions: Dict[Any, FlagResolution] = {}
        self.flags: Dict[Any, Flag] = {}
        self.flag_comments: Dict[Any, CommentDetails] = {}
        self.users: Dict[int, User] = {}
        self.created_flag_id: Optional[int] = None

        self._include_collections = False
        self._include_flag_types = False

        contexts = set(collection.context for collection in self.flag_collections)
        for flag_type in FlagType.objects.filter(context__in=contexts):
            self.include_flag_type(flag_type)

        for resolution in FlagResolution.objects.all():
            self.flag_resolutions[resolution.id] = resolution

    def include_just_closed(self):
        now = timezone.now()
        yesterday = now - datetime.timedelta(minutes=2)
        recently_closed_ids = FlagComment.objects.filter(
            created__gte=yesterday,
            flag__collection__in=self.flag_collections).exclude(flag__resolution__status=FlagStatus.OPEN)\
            .values_list('flag', flat=True).distinct()

        recently_cosed_flags = Flag.objects.filter(pk__in=recently_closed_ids).select_related('user', 'collection')
        for flag in recently_cosed_flags:
            self.include_flag(flag)

    def include_user(self, user: User):
        self.users[user.id] = user

    def add_flag(self, data: dict, flag_collection: FlagCollection = None):
        if not flag_collection:
            flag_collection = self.flag_collections[0]

        flag_type = data.pop('flag_type')
        comment = data.pop('comment', None)
        user_private = data.pop('user_private', False)
        resolution = data.pop('resolution', None)
        if resolution:
            resolution = FlagResolution.objects.get(pk=resolution)

        flag = flag_collection.add_flag(
            FlagType.objects.get(pk=flag_type),
            user=self.user,
            comment=comment,
            user_private=user_private,
            resolution=resolution
        )
        self.include_flag(flag)
        self.created_flag_id = flag.id
        return self

    def include_flag_types(self):
        self._include_flag_types = True
        return self

    def include_open_flags(self):
        # default is to get open flags
        #for flag_collection in self.flag_collections:
        #    for flag in flag_collection.flags(user=self.user).filter(FlagCollection.Q_OPEN_FLAGS):
        #        self.include_flag(flag)
        flags = Flag.objects.filter(collection__in=self.flag_collections).filter(FlagCollection.Q_OPEN_FLAGS)\
            .select_related('user', 'collection')
        for flag in flags:
            self.include_flag(flag)

        return self

    def include_comment(self, flag_comment: FlagComment, detailed: bool = False):
        if flag_comment.id not in self.flag_comments:
            self.flag_comments[flag_comment.id] = CommentDetails(helper=self, comment=flag_comment)
        fc = self.flag_comments[flag_comment.id]
        fc.detailed = fc.detailed or detailed
        self.include_user(flag_comment.user)
        if flag_comment.flag_id not in self.flags:
            self.include_flag(flag_comment.flag, include_comments=False)

    def include_comments_since(self, since):
        comments = FlagComment.objects.filter(flag__collection__in=self.flag_collections, created__gt=since).order_by('-created').select_related('flag', 'flag__collection')
        comments = [comment for comment in comments if self.is_viewable_flag(comment.flag)]
        for comment in comments:
            if comment.created > since:
                since = comment.created

            if self.is_viewable_flag(comment.flag):
                self.include_comment(comment, detailed=True)
        self.since = since

    def is_viewable_flag(self, flag: Flag) -> bool:
        if flag.user_private and \
            flag.user != self.user and \
            not flag.collection.is_owner_or_admin(self.user):

            return False
        return True

    def include_history(self, collection):
        flags = Flag.objects.filter(collection=collection)
        for flag in flags:
            self.include_flag(flag, include_comments=True)

    def include_flag(self, flag: Flag, include_comments: bool = False):

        if not self.is_viewable_flag(flag):
            #print('debug, filtering out private flag')
            return

        self.flags[flag.id] = flag

        if include_comments:
            comments = flag.flagcomment_set.order_by('created')
            for comment in comments:
                self.include_comment(flag_comment=comment, detailed=True)
            # previously included first and last comment for tooltips
            # but made syncing difficult
            """
            else:
                if comments:
                    comments = comments.all()
                    self.include_comment(flag_comment=comments.first())
                    self.include_comment(flag_comment=comments.last())
            """
        self.include_user(flag.user)

    def include_flag_type(self, flag_type: FlagType):
        self.flag_types[flag_type.id] = flag_type

    def include_collections(self):
        self._include_collections = True

    def lab_text(self, user: User) -> str:
        if user.is_superuser:
            return 'admin'
        labs = list(Lab.valid_labs_qs(user))
        if len(labs) == 0:
            return 'no affiliation'
        if len(labs) == 1:
            return labs[0].name

        orgs = set()
        for lab in labs:
            orgs.add(lab.organization.name)
        orgs = list(orgs)
        orgs.sort()
        return ', '.join(orgs)

    def to_json(self):
        json_data = {}
        flags_json = []

        users_json = []
        for user in self.users.values():
            user_settings = UserSettings.get_for_user(user)
            json_entry = {
                'id': user.id,
                'name': UserSettings.preferred_label_for(user),
                'avatar': user_settings.avatar_url,
                'color': user_settings.avatar_color,
                'lab': self.lab_text(user)
            }
            users_json.append(json_entry)
        json_data['users'] = users_json

        flag_collections_json = []
        flags_info = None
        if self._include_collections:

            flags_info = fetch_flag_infos(flag_collections=self.flag_collections, flags=self.flags.values(), user=self.user)

            for flag_collection in self.flag_collections:
                permission_level = flag_collection.permission_level(self.user)

                json_entry = {
                    'id': flag_collection.id,
                    'user_permission': permission_level.level,
                    'context': flag_collection.context_id,
                    # TODO watching isn't currently used
                    #'watching': flag_collection.unseen_flag_activity(user=self.user)
                }
                json_entry = {**json_entry, **flag_collection.extra_info}
                flag_collections_json.append(json_entry)

        json_data['collections'] = flag_collections_json

        for flag in self.flags.values():
            flag_type = self.flag_types[flag.flag_type_id]
            if flag_type.only_one:
                # find when flag was last opened
                try:
                    created = FlagComment.objects.filter(flag=flag, resolution__status=FlagStatus.OPEN).order_by('-created').\
                        values_list('created', flat=True).first()
                except:
                    pass

            resolution = self.flag_resolutions[flag.resolution_id]
            json_entry = {
                'id': flag.pk,
                'collection': flag.collection_id,
                'flag_type': flag.flag_type_id,
                'status': resolution.status,
                'resolution': resolution.id,
                'user': flag.user_id,
                'created': flag.created.timestamp()
            }
            if flag.data:
                json_entry['data'] = flag.data
            if flags_info:
                extra_flag_data = flags_info.extra_flag_info(flag)
                # want this set toa  dict or None (None so we override any previous extra data)
                json_entry['extra_data'] = extra_flag_data

            flags_json.append(json_entry)
        json_data['flags'] = flags_json

        if self._include_flag_types:
            flag_types = []

            for flag_type in self.flag_types.values():
                resolutions = list(FlagTypeResolution.objects.filter(flag_type=flag_type).values_list('resolution__id', flat=True))
                json_entry = {
                    'id': flag_type.pk,
                    'only_one': flag_type.only_one,
                    'context': flag_type.context_id,
                    'label': flag_type.label,
                    'description': flag_type.description,
                    'help_text': flag_type.help_text,
                    'permission': flag_type.permission_enum.level,
                    'raise_permission': flag_type.raise_permission_enum.level,
                    'comments_enabled': flag_type.comments_enabled,
                    'importance': flag_type.importance,
                    'resolutions': resolutions
                }

                flag_types.append(json_entry)
            json_data['flag_types'] = flag_types

            flag_resolutions = []
            for flag_resolution in self.flag_resolutions.values():
                json_entry = {
                    'id': flag_resolution.id,
                    'label': flag_resolution.label,
                    'status': flag_resolution.status
                }
                flag_resolutions.append(json_entry)
            json_data['flag_resolutions'] = flag_resolutions
        #end include flag types

        flag_comments = []
        for comment in self.flag_comments.values():
            json_entry = comment.to_json()
            flag_comments.append(json_entry)
        json_data['comments'] = flag_comments
        if self.since:
            json_data['since'] = self.since.timestamp()

        if self.created_flag_id:
            json_data['created_flag_id'] = self.created_flag_id

        return json_data


class FlagsView(APIView):

    def get(self, request, **kwargs) -> Response:
        pks = kwargs.get('flag_collection_id').split(',')
        fcs = FlagCollection.objects.filter(pk__in=pks)
        fcs = fcs.select_related("context")

        flag_helper = FlagHelper(flag_collections=fcs, user=request.user)

        # history already assumes you got the non-history part with collections and flag types
        history = request.query_params.get('history')
        since = request.query_params.get('since')
        flag_helper.since = ensure_timezone_aware(datetime.datetime.utcnow())

        if history:
            # user has now seen the flags, mark them as such
            history = int(history)
            FlagWatch.objects.filter(flag_collection=history, user=request.user).update(since=now())
            flag_helper.include_history(history)

        if since:
            since = ensure_timezone_aware(datetime.datetime.fromtimestamp(float(since)))
            flag_helper.include_comments_since(since)

        if not history and not since:
            flag_helper.include_flag_types()
            flag_helper.include_collections()
            flag_helper.include_open_flags()
            flag_helper.include_just_closed()

        flag_helper.include_user(request.user)

        return Response(flag_helper.to_json())

    def post(self, request, **kwargs) -> Response:
        pk = kwargs.get('flag_collection_id')
        fc = FlagCollection.objects.get(pk=pk)

        flag_helper = FlagHelper(flag_collections=fc, user=request.user)

        watch = request.data.get('watch')
        if watch is not None:
            fc.set_watcher(user=request.user, watch=watch)
            flag_helper.include_collections()
        else:
            flag_helper.add_flag(data=request.data)

        return Response(flag_helper.to_json())


class FlagView(APIView):

    def get(self, request, **kwargs) -> Response:
        pk = kwargs.get('flag_id')
        f = Flag.objects.get(pk=pk)
        # TODO - add comments
        flag_helper = FlagHelper(flag_collections=f.collection, user=request.user)
        flag_helper.include_flag(f, include_comments=True)
        flag_helper.include_user(request.user)
        return Response(flag_helper.to_json())

    def post(self, request, **kwargs) -> Response:
        pk = kwargs.get('flag_id')
        f = Flag.objects.get(pk=pk)  # type: Flag
        data = request.data
        resolution = data.get('resolution', None)
        if resolution:
            resolution = FlagResolution.objects.get(pk=resolution)

        f.flag_action(
            user=request.user,
            resolution=resolution,
            comment=empty_to_none(data.get('comment')))

        flag_helper = FlagHelper(flag_collections=f.collection, user=request.user)
        flag_helper.include_flag(f, include_comments=True)
        return Response(flag_helper.to_json())
