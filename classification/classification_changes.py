from datetime import datetime
from typing import Any, Iterator, List, Optional, Iterable, Union

from django.db.models import QuerySet, Q
from django.contrib.auth.models import User
from flags.models import FlagComment, Flag, FlagResolution
from snpdb.models import Allele
from classification.enums import SubmissionSource
from classification.models import ClassificationModification, Classification, \
    classification_flag_types
from django.utils.timezone import now


class InfiniteQuerySubject:

    @staticmethod
    def chunk_query(qs: QuerySet, chunk_size=50) -> Iterable[Any]:
        skip = 0
        while True:
            cache = list(qs.all()[skip: skip+chunk_size])
            skip = skip + chunk_size
            if not cache:
                return None
            for val in cache:
                yield val

    def __init__(self, qs: QuerySet):
        self.iterator = InfiniteQuerySubject.chunk_query(qs)
        self.subject = None
        self.finished = False
        self.next()

    def next(self):
        while True:
            if self.finished:
                return
            subject = None
            try:
                subject = next(self.iterator)
            except StopIteration:
                self.finished = True
                return

            changes = ClassificationChanges.from_object(subject)
            if changes and changes.changes:
                self.subject = changes
                return


class Stitch:
    def __init__(self, queries: List[QuerySet]):
        self.subjects = [InfiniteQuerySubject(qs) for qs in queries]

    def iterate(self) -> Iterator['ClassificationChanges']:
        while True:
            min_iqs = None
            for iqs in self.subjects:
                if not iqs.finished:
                    if min_iqs is None:
                        min_iqs = iqs
                    elif iqs.subject.date > min_iqs.subject.date:
                        min_iqs = iqs
            if not min_iqs:
                return None
            subject = min_iqs.subject
            min_iqs.next()
            yield subject


class ClassificationChange:
    def __init__(self, key: str, attribute: str, before: Any, after: Any, resolution: Optional[FlagResolution] = None):
        self.key = key
        self.attribute = attribute
        self.before = before
        self.after = after
        self.resolution = resolution

    def __str__(self):
        return f'{self.key}.{self.attribute}'


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
    def list_changes(classification: Optional[Classification] = None, latest_date: Optional[datetime] = None, limit=200) -> List['ClassificationChanges']:
        if not latest_date:
            latest_date = now()

        date_q = Q(created__lte=latest_date)

        vcm_q = date_q
        flags_q = date_q

        if classification:
            vcm_q = vcm_q & Q(classification=classification)
            flags_q = flags_q & Q(flag__collection=classification.flag_collection_safe)
        else:
            pass

        # Variant Classification Changes
        vcm_qs = ClassificationModification.objects.filter(vcm_q) \
                     .select_related('classification', 'classification__lab',
                                     'classification__user').order_by('-created')

        # Flag Changes
        flags_qs = FlagComment.objects.filter(
            flags_q
        ).exclude(flag__flag_type__in=[
            classification_flag_types.classification_outstanding_edits,
            classification_flag_types.unshared_flag
        ]).select_related('flag', 'flag__flag_type', 'flag__collection', 'resolution').order_by('-created')

        stitch = Stitch([vcm_qs, flags_qs])
        changes = []
        for index, vcmc in enumerate(stitch.iterate()):
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

    def __str__(self):
        return f'{self.date} - {" ".join([str(c) for c in self.changes])}'
