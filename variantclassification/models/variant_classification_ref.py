from datetime import datetime, timezone
from typing import Optional

from django.contrib.auth.models import User
from django.core.exceptions import PermissionDenied
from django.http.response import Http404
from lazy import lazy
import re

from library.utils import empty_to_none
from snpdb.models import Lab
from variantclassification.enums import SubmissionSource
from variantclassification.models.variant_classification import VariantClassification, \
    VariantClassificationProcessError, VariantClassificationModification
from variantclassification.models import VariantClassificationJsonParams


class VariantClassificationRef:

    ALLOWED_LAB_ID_RE = re.compile('^[A-Za-z0-9_-]*$')

    """
    Typically used to change a string reference of a variant classification record to match the actual record
    """

    def __init__(self, user: User, rid: int = None, lab: Lab = None, lab_record_id: str = None, version: int = None):
        """
        @param rid the record id (the normal one)
        """
        self.user = user
        self.rid = rid
        self.lab = lab
        self.lab_record_id = lab_record_id
        self.version = version
        self.cached_version = None
        self.cached_record = None

    @lazy
    def record(self) -> VariantClassification:
        """
        if init_from_obj we'll have a cached record
        """
        if self.cached_record:
            return self.cached_record

        if self.rid:
            return VariantClassification.objects.filter(pk=self.rid).first()

        if self.lab_record_id:
            return VariantClassification.objects.filter(lab=self.lab, lab_record_id=self.lab_record_id).first()

    @lazy
    def modification(self) -> VariantClassificationModification:
        _record = self.record
        if self.version_datetime:
            vcm = VariantClassificationModification.objects\
                .select_related('variant_classification', 'variant_classification__lab')\
                .get(variant_classification=self.record.id, created=self.version_datetime)
            vcm.check_can_view(self.user)
            return vcm
        return VariantClassificationModification.latest_for_user(self.user, self.record, exclude_withdrawn=False).first()

    def check_exists(self):
        if not self.exists():
            raise Http404(f"Could not find record {self.rid}")

    def check_security(self, must_be_writable=False):
        if self.record is None and not must_be_writable:
            raise VariantClassificationProcessError('Version ' + str(self.version) + ' not found')
        if self.record is not None:
            if self.version is None:

                # use can write to the record, so no need to version the record
                if self.record.can_write(self.user):
                    return
                if must_be_writable:
                    raise PermissionDenied()

                latest_version = self.record.latest_modification_for_user(self.user, exclude_withdrawn=False)
                if not latest_version:
                    raise PermissionDenied()
                self.version = latest_version.created.timestamp()
            else:
                if must_be_writable:
                    raise VariantClassificationProcessError('Can not perform actions on a version')

                mod = VariantClassificationModification.objects.filter(variant_classification=self.record, created=self.version_datetime).first()

                if mod is None:
                    raise VariantClassificationProcessError('Version ' + str(self.version) + ' not found')
                if not mod.can_view(self.user):
                    raise PermissionDenied()

    @lazy
    def version_datetime(self) -> Optional[datetime]:
        if self.version:
            return datetime.utcfromtimestamp(self.version).replace(tzinfo=timezone.utc)
        return None

    def exists(self) -> bool:
        return self.record is not None

    def create(self,
               source: SubmissionSource,
               data: dict = None,
               save: bool = True,
               make_fields_immutable=False):
        if not self.lab:
            raise ValueError('Cannot create record without a lab')
        if self.exists():
            raise VariantClassificationProcessError("Can't create a record when one already exists")
        if self.rid:
            raise VariantClassificationProcessError("Can't create a record from a plain ID, must specific lab and lab record id")
        if not VariantClassificationRef.ALLOWED_LAB_ID_RE.match(self.lab_record_id):
            raise VariantClassificationProcessError(f"lab_id \"{self.lab_record_id}\" contains illegal characters - only letters, numbers, underscores and dashes allowed")

        self.cached_record = VariantClassification.create(
            user=self.user,
            lab=self.lab,
            lab_record_id=self.lab_record_id,
            data=data,
            save=save,
            source=source,
            make_fields_immutable=make_fields_immutable,
        )
        return self.cached_record

    def as_json(self,
                params: VariantClassificationJsonParams) -> dict:

        params.version = self.cached_version or self.version
        params.include_data = params.include_data or not params.version is None

        return self.record.as_json(params)

    @staticmethod
    def init_from_obj(user: User, obj) -> 'VariantClassificationRef':
        if isinstance(obj, VariantClassificationRef):
            return obj
        if isinstance(obj, VariantClassification):
            vc = VariantClassificationRef(user=user)
            vc.cached_record = obj
            return vc
        if isinstance(obj, VariantClassificationModification):
            vc = VariantClassificationRef(user=user)
            vc.cached_record = obj.variant_classification
            vc.cached_version = obj
            vc.version = obj.created.timestamp()
            return vc
        if isinstance(obj, str):
            return VariantClassificationRef.init_from_str(user, obj)
        if isinstance(obj, dict):
            rid = obj.get('id')
            return VariantClassificationRef.init_from_obj(user, VariantClassificationModification.objects.get(pk=rid))

        raise VariantClassificationProcessError("Can not convert " + str(obj) + " to a VariantClassificationRef")

    @staticmethod
    def init_from_parts(user: User, lab_id: str = None, record_id: str = None, lab_record_id: str = None, version: str = None) -> 'VariantClassificationRef':
        if version:
            version = float(version)

        # don't have access to this group
        lab = None
        if lab_id:
            if not (user.is_superuser or user.groups.filter(name=lab_id).exists()):
                raise VariantClassificationProcessError('Submitting user does not have access to lab ' + lab_id)
            lab = Lab.objects.filter(group_name=lab_id).first()
            if not lab:
                raise VariantClassificationProcessError(f"No such lab as ({lab_id})")
        elif lab_record_id:
            raise VariantClassificationProcessError('Can provide lab_record_id without a lab')

        return VariantClassificationRef(
            user=user, lab=lab, rid=record_id, lab_record_id=lab_record_id, version=version)

    @staticmethod
    def init_from_str(user: User, id_str: str) -> Optional['VariantClassificationRef']:
        if not id_str:
            return None
        id_str = str(id_str)

        # group 1 is the lab
        # group 2 is the lab_record_id or regular id
        # group 3 is the version
        match = re.match(r'^(?:([^.]*?)(?:/)?)([^/.]*?)(?:[.](.*))?$', id_str)

        lab_ref = empty_to_none(match[1])
        record_id = empty_to_none(match[2])
        version = empty_to_none(match[3])

        kwargs = {}
        if lab_ref:
            kwargs["lab_id"] = lab_ref
            kwargs["lab_record_id"] = record_id
        else:
            kwargs["record_id"] = record_id

        return VariantClassificationRef.init_from_parts(user=user, version=version, **kwargs)
