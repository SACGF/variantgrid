import operator
from datetime import datetime
from functools import total_ordering
from typing import Any, List, Optional, Union

from django.contrib.auth.models import User
from django.db.models import Q, Subquery
from django.utils.timezone import now

from classification.enums import SubmissionSource
from classification.models import ClassificationModification, Classification, classification_flag_types
from flags.models import FlagComment, Flag, FlagResolution
from library.utils import IterableTransformer, IteratableStitcher
from snpdb.models import Allele, Lab


class ClassificationChange:
    def __init__(self, key: str, attribute: str, before: Any, after: Any, resolution: Optional[FlagResolution] = None):
        self.key = key
        self.attribute = attribute
        self.before = before
        self.after = after
        self.resolution = resolution

    def __str__(self):
        return f'{self.key}.{self.attribute}'

    @property
    def is_large(self):
        if isinstance(self.before, str) and (len(self.before) > 50 or '\n' in self.before):
            return True
        if isinstance(self.after, str) and (len(self.after) > 50 or '\n' in self.after):
            return True
        return False


@total_ordering
class ClassificationChanges:

    @staticmethod
    def from_object(obj) -> Optional['ClassificationChanges']:
        if isinstance(obj, FlagComment):
            return ClassificationChanges.from_flag_comment(obj)
        if isinstance(obj, ClassificationModification):
            return ClassificationChanges.from_vcm(obj)
        raise ValueError(f'Cannot convert {obj} to ClassificationChanges')

    @staticmethod
    def from_flag_comment(fc: FlagComment) -> Optional['ClassificationChanges']:
        flag = fc.flag
        collection = flag.collection
        record = Classification.objects.filter(flag_collection=collection).first()
        if not record:
            record = Allele.objects.filter(flag_collection=collection).first()

        if not record:
            return None

        return ClassificationChanges(
            record=record,
            date=fc.created,
            user=fc.user,
            changes=[ClassificationChange(
                key='Flag',
                attribute=flag.flag_type.label,
                before=None,
                after=fc.text,
                resolution=fc.resolution
            )],
            is_creation=False,  # TODO
            ignore_before=True,
            flag=flag)

    @staticmethod
    def from_vcm(vcm: ClassificationModification, include_debug=False) -> 'ClassificationChanges':
        prev = vcm.previous
        is_creation = prev is None

        evidence = prev.evidence if prev else {}
        delta = vcm.delta

        changes: List[ClassificationChange] = []

        def append_change(key, attribute, original_val, value):
            if original_val == value:
                return
            nonlocal include_debug
            nonlocal changes
            if include_debug or attribute in ('value', 'note', 'explain'):
                changes.append(ClassificationChange(key, attribute, original_val, value))

        for key, blob in delta.items():
            original = evidence.get(key, {})
            if blob:
                for attribute, value in blob.items():
                    original_val = original.get(attribute)
                    append_change(key, attribute, original_val, value)
            else:
                for attribute, original_val in original.items():
                    append_change(key, attribute, original_val, None)

        changes.sort(key=lambda c: f'{c.key}-{c.attribute}')

        return ClassificationChanges(
            record=vcm.classification,
            date=vcm.created,
            user=vcm.user,
            changes=changes,
            is_creation=is_creation,
            ignore_before=is_creation,
            source=vcm.source)

    @staticmethod
    def list_changes(
            classifications: Optional[List[Classification]] = None,
            latest_date: Optional[datetime] = None,
            for_user: Optional[User] = None,
            for_lab: Optional[Lab] = None,
            limit=100) -> List['ClassificationChanges']:
        if not latest_date:
            latest_date = now()

        date_q = Q(created__lte=latest_date)

        vcm_q = date_q
        flags_q = date_q

        if classifications:
            cids = [c.pk for c in classifications]
            flag_collections = [fid for fid in [c.flag_collection_id for c in classifications] if fid]
            vcm_q = vcm_q & Q(classification__in=cids)
            flags_q = flags_q & Q(flag__collection__in=flag_collections)
        elif for_user:
            vcm_q = vcm_q & Q(user=for_user)
            flags_q = flags_q & Q(user=for_user)
        elif for_lab:
            vcm_q = vcm_q & Q(classification__lab=for_lab)

            flag_collections_q = Classification.objects.filter(lab=for_lab).values_list('flag_collection', flat=True)
            flags_q = flags_q & Q(flag__collection__context__id='classification')
            flags_q = flags_q & Q(flag__collection__in=Subquery(flag_collections_q))
        else:
            # can't currently render clinical context flags
            flags_q = flags_q & Q(flag__collection__context__id__in=['classification', 'allele'])

        # Variant Classification Changes
        vcm_qs = ClassificationModification.objects.filter(vcm_q) \
                     .select_related('classification', 'classification__lab',
                                     'classification__user', 'user').order_by('-created')[:limit+1]

        # Flag Changes
        flags_qs = FlagComment.objects.filter(
            flags_q
        ).exclude(flag__flag_type__in=[
            classification_flag_types.classification_outstanding_edits,
            classification_flag_types.unshared_flag
        ]).select_related('flag', 'flag__flag_type', 'flag__collection', 'resolution', 'user').order_by('-created')[:limit+1]

        stitcher = IteratableStitcher(
            iterables=[
                IterableTransformer(vcm_qs, ClassificationChanges.from_object),
                IterableTransformer(flags_qs, ClassificationChanges.from_object)
            ],
            comparison=operator.__gt__
        )

        changes = []
        for index, vcmc in enumerate(stitcher):
            changes.append(vcmc)
            if index >= limit:
                break

        return changes

    def __init__(self,
                 record: Union[Classification, Allele],
                 date: datetime,
                 user: User,
                 changes: List[ClassificationChange],
                 is_creation: bool,
                 ignore_before: bool = False,
                 source: str = None,
                 flag: Flag = None):
        self.record = record
        self.date = date
        self.user = user
        self.changes = changes
        self.ignore_before = ignore_before
        self.is_creation = is_creation
        self.source = source
        self.flag = flag

    @property
    def classification(self) -> Optional[Classification]:
        if isinstance(self.record, Classification):
            return self.record
        return None

    @property
    def allele(self) -> Optional[Allele]:
        if isinstance(self.record, Allele):
            return self.record
        return None

    @property
    def record_id(self):
        vc = self.classification
        if vc:
            return f'vc{vc.id}'

        allele = self.allele
        if allele:
            return f'a{allele.id}'
        return None

    @property
    def source_label(self) -> str:
        if self.flag:
            return self.flag.flag_type.label
        source_label = {SubmissionSource.API: 'API',
                        SubmissionSource.VARIANT_GRID: 'Autopopulate',
                        SubmissionSource.CONSENSUS: 'Copy from latest classification',
                        SubmissionSource.FORM: 'Form entry'}
        return source_label.get(self.source, self.source)

    def __lt__(self, other):
        return self.date < other.date

    def __str__(self):
        return f'{self.date} - {" ".join([str(c) for c in self.changes])}'
