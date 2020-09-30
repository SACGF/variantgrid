from datetime import datetime
from typing import Any, Iterator, List, Optional, Iterable, Union

from django.db.models import QuerySet, Q
from django.contrib.auth.models import User
from flags.models import FlagComment, Flag, FlagResolution
from snpdb.models import Allele
from variantclassification.enums import SubmissionSource
from variantclassification.models import VariantClassificationModification, VariantClassification, \
    variant_classification_flag_types
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

            changes = VariantClassificationChanges.from_object(subject)
            if changes and changes.changes:
                self.subject = changes
                return


class Stitch:
    def __init__(self, queries: List[QuerySet]):
        self.subjects = [InfiniteQuerySubject(qs) for qs in queries]

    def iterate(self) -> Iterator['VariantClassificationChanges']:
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


class VariantClassificationChange:
    def __init__(self, key: str, attribute: str, before: Any, after: Any, resolution: Optional[FlagResolution] = None):
        self.key = key
        self.attribute = attribute
        self.before = before
        self.after = after
        self.resolution = resolution

    def __str__(self):
        return f'{self.key}.{self.attribute}'


class VariantClassificationChanges:

    @staticmethod
    def from_object(obj) -> Optional['VariantClassificationChanges']:
        if isinstance(obj, FlagComment):
            return VariantClassificationChanges.from_flag_comment(obj)
        if isinstance(obj, VariantClassificationModification):
            return VariantClassificationChanges.from_vcm(obj)
        raise ValueError(f'Cannot convert {obj} to VariantClassificationChanges')

    @staticmethod
    def from_flag_comment(fc: FlagComment) -> Optional['VariantClassificationChanges']:
        flag = fc.flag
        collection = flag.collection
        record = VariantClassification.objects.filter(flag_collection=collection).first()
        if not record:
            record = Allele.objects.filter(flag_collection=collection).first()

        if not record:
            return None

        return VariantClassificationChanges(
            record=record,
            date=fc.created,
            user=fc.user,
            changes=[VariantClassificationChange(
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
    def from_vcm(vcm: VariantClassificationModification, include_debug=False) -> 'VariantClassificationChanges':
        prev = vcm.previous
        is_creation = prev is None

        evidence = prev.evidence if prev else {}
        delta = vcm.delta

        changes: List[VariantClassificationChange] = []

        def append_change(key, attribute, original_val, value):
            if original_val == value:
                return
            nonlocal include_debug
            nonlocal changes
            if include_debug or attribute in ('value', 'note', 'explain'):
                changes.append(VariantClassificationChange(key, attribute, original_val, value))

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

        return VariantClassificationChanges(
            record=vcm.variant_classification,
            date=vcm.created,
            user=vcm.user,
            changes=changes,
            is_creation=is_creation,
            ignore_before=is_creation,
            source=vcm.source)

    @staticmethod
    def list_changes(variant_classification: Optional[VariantClassification] = None, latest_date: Optional[datetime] = None, limit=200) -> List['VariantClassificationChanges']:
        if not latest_date:
            latest_date = now()

        date_q = Q(created__lte=latest_date)

        vcm_q = date_q
        flags_q = date_q

        if variant_classification:
            vcm_q = vcm_q & Q(variant_classification=variant_classification)
            flags_q = flags_q & Q(flag__collection=variant_classification.flag_collection_safe)
        else:
            pass

        # Variant Classification Changes
        vcm_qs = VariantClassificationModification.objects.filter(vcm_q) \
                     .select_related('variant_classification', 'variant_classification__lab',
                                     'variant_classification__user').order_by('-created')

        # Flag Changes
        flags_qs = FlagComment.objects.filter(
            flags_q
        ).exclude(flag__flag_type__in=[
            variant_classification_flag_types.variant_classification_outstanding_edits,
            variant_classification_flag_types.unshared_flag
        ]).select_related('flag', 'flag__flag_type', 'flag__collection', 'resolution').order_by('-created')

        stitch = Stitch([vcm_qs, flags_qs])
        changes = []
        for index, vcmc in enumerate(stitch.iterate()):
            changes.append(vcmc)
            if index >= limit:
                break

        return changes

    def __init__(self,
                 record: Union[VariantClassification, Allele],
                 date: datetime,
                 user: User,
                 changes: List[VariantClassificationChange],
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
    def variant_classification(self) -> Optional[VariantClassification]:
        if isinstance(self.record, VariantClassification):
            return self.record
        return None

    @property
    def allele(self) -> Optional[Allele]:
        if isinstance(self.record, Allele):
            return self.record
        return None

    @property
    def record_id(self):
        vc = self.variant_classification
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
