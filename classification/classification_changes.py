import operator
from datetime import datetime
from functools import total_ordering
from typing import Any, Optional, Union

from django.contrib.auth.models import User
from django.db.models import Q, Subquery
from django.utils.timezone import now

from classification.enums import SubmissionSource
from classification.models import ClassificationModification, Classification, classification_flag_types
from flags.models import FlagComment, Flag, FlagResolution
from library.utils import IterableTransformer, IterableStitcher
from snpdb.models import Allele, Lab


class ClassificationChange:
    """
    For showing on Classification Activity Page.
    e.g. "PM1" changed "PM" -> "PS"
    ClassificationChange will be bundled into ClassificationChanges for a classification ID, date and series of changes
    Changes include newly created records
    """

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
    """
    Captures all changes to a Classification from one action (upload, web form edits, creation from web form or flag changing).
    May involve 1 to many ClassificationChange (each change representing a different EvidenceKey if multiple EvidenceKeys were chagned at once)
    """

    def __init__(self,
                 record: Union[Classification, Allele],
                 date: datetime,
                 user: User,
                 changes: list[ClassificationChange],
                 is_creation: bool,
                 ignore_before: bool = False,
                 source: str = None,
                 flag: Flag = None):
        """
        :param record: The Classification that has been modified (TODO do we ever pass Alleles to this? It would just be for Flags affecting Alleles which is not actively used)
        :param date: The datetime the change took place
        :param user: The user that initiated the change
        :param changes: A full list of changes that happened during this change
        :param is_creation: Is the "change" actually the creation of the record
        :param ignore_before: Instead of showing a value going from X->Y, just say it is Y (typically should be True if is_creation or flag). TODO auto calculate based on those other params.
        :param source: The source, typically a a SubmissionSource
        :param flag: If a flag is the cause of the change, provide it here
        """
        self.record = record
        self.date = date
        self.user = user
        self.changes = changes
        self.ignore_before = ignore_before
        self.is_creation = is_creation
        self.source = source
        self.flag = flag

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

        changes: list[ClassificationChange] = []

        def append_change(key: str, attribute: str, original_val: Any, value: Any) -> None:
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
            for_classifications: Optional[list[Classification]] = None,
            for_user: Optional[User] = None,
            for_lab: Optional[Lab] = None,
            latest_date: Optional[datetime] = None,
            limit=100) -> list['ClassificationChanges']:
        """
        Create a list of ClassificationChanges filtered by the given parameters.
        Will filter on classifications (if not provided then) filter for_user (if not provided then) filter for for_lab.
        NOTE filters don't combine

        :param for_classifications: An optional list of classifications to only look at.
        :param for_user: If provided, only changes performed by the user
        :param for_lab: If provided only changes against Classifications for the given lab
        :param latest_date: Get limit changes from latest_date or before (handy for showing most recent changes).
        :param limit: The max number of changes to list
        :return: A list of ClassificationChanges ordered latest_date to earlier, no longer than limit long
        """

        if not latest_date:
            latest_date = now()

        date_q = Q(created__lte=latest_date)

        vcm_q = date_q
        flags_q = date_q

        if for_classifications:
            c_ids = [c.pk for c in for_classifications]
            flag_collections = [fid for fid in [c.flag_collection_id for c in for_classifications] if fid]
            vcm_q &= Q(classification__in=c_ids)
            flags_q &= Q(flag__collection__in=flag_collections)
        elif for_user:
            vcm_q &= Q(user=for_user)
            flags_q &= Q(user=for_user)
        elif for_lab:
            vcm_q &= Q(classification__lab=for_lab)

            flag_collections_q = Classification.objects.filter(lab=for_lab).values_list('flag_collection', flat=True)
            flags_q &= Q(flag__collection__context__id='classification')
            flags_q &= Q(flag__collection__in=Subquery(flag_collections_q))
        else:
            # can't currently render clinical context flags
            flags_q &= Q(flag__collection__context__id__in=['classification', 'allele'])

        # Classification Changes
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

        stitcher = IterableStitcher(
            iterables=[
                IterableTransformer(vcm_qs, ClassificationChanges.from_object),
                IterableTransformer(flags_qs, ClassificationChanges.from_object)
            ],
            comparison=operator.__gt__
        )

        # both classification changes and flag changes have been limited to "limit", but when combined there could be "limit x 2" records
        # so reduce down to "limit"
        changes = []
        for index, vcmc in enumerate(stitcher):
            changes.append(vcmc)
            if index >= limit:
                break

        return changes

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
