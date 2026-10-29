"""
Patients and the material taken from them: Patient (phenotype text matched to ontology terms, which
patients/models/models_phenotype.py holds),
Specimen (one tissue at one timepoint), Extraction (nucleic acid off a specimen) and
SpecimenMeasure, all Guardian-permissioned and optionally externally managed (ExternalPK /
ExternallyManagedModel). ExtractionMatchMixin is how a Sample claims its extraction before the
records exist. Patient modifications and imports are audited rows; clinicians and patient records
complete the set.
"""
import datetime
import operator
import os
from functools import reduce
from typing import Optional

import nameparser
from django.conf import settings
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.db import models, transaction
from django.db.models import Case, F, Q, QuerySet, TextField, Value, When
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.functions import Concat
from django.dispatch.dispatcher import Signal, receiver
from django.urls.base import reverse
from django.utils import timezone
from easy_thumbnails.files import get_thumbnailer
from guardian.models import GroupObjectPermission, UserObjectPermission
from model_utils.models import TimeStampedModel

from library.django_utils import ensure_mutally_exclusive_fields_not_set
from library.django_utils.django_file_system_storage import PrivateUploadStorage
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.enums.file_attachments import AttachmentFileType
from library.enums.titles import Title
from library.preview_request import PreviewData, PreviewKeyValue, PreviewModelMixin
from library.utils import calculate_age
from patients.enums import (
    MatchStatus,
    NucleicAcid,
    PatientRecordMatchType,
    PatientRecordOriginType,
    PopulationGroup,
    Sex,
    SpecimenMeasureType,
    TissueStatus,
)
from patients.external_references import ResolvedReference, normalize_specimen_reference
from patients.models.has_phenotype_description_mixin import HasPhenotypeDescriptionMixin

TEST_PATIENT_KWARGS = {"first_name": "PATIENT", "last_name": "TESTPATIENT"}

# Patient to the ontology terms matched in its phenotype text
PATIENT_ONTOLOGY_TERM_PATH = "phenotype_description__textphenotypesentence__text_phenotype__textphenotypematch__ontology_term"


class FakeData(TimeStampedModel):
    pass


class ExternalModelManager(TimeStampedModel):
    """ A way to say a model is externally managed """
    name = models.TextField(primary_key=True)
    details = models.TextField(blank=True)
    can_modify = models.BooleanField(default=False)  # People in VG can modify the model (via forms etc)
    modifications_sent_to_external_system = models.BooleanField(default=False)

    @property
    def explaination(self):
        msg = f"This object is a copy of data from the external system {self}."
        if not self.can_modify:
            msg += " You cannot modify it."
        else:
            if self.modifications_sent_to_external_system:
                msg += " Modifications will be sent to external system."
        return msg

    def __str__(self):
        name = self.name
        if self.details:
            name += f" ({self.details})"
        return name


class ExternalPK(TimeStampedModel, PreviewModelMixin):
    code = models.TextField()
    external_type = models.TextField()
    external_manager = models.ForeignKey(ExternalModelManager, on_delete=CASCADE)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-person-walking-arrow-right"

    @classmethod
    def preview_if_url_visible(cls) -> Optional[str]:
        return 'patients'

    class Meta:
        unique_together = ('code', 'external_type', 'external_manager')

    def __str__(self):
        return f"{self.code} ({self.external_type})"


class ExternallyManagedModel(TimeStampedModel):
    external_pk = models.OneToOneField(ExternalPK, null=True, on_delete=CASCADE)
    # The field holding the local reference beside external_pk - @see patients.external_references
    LOCAL_REFERENCE_FIELD = "reference_id"

    class Meta:
        abstract = True

    @classmethod
    def normalize_reference_id(cls, reference_id: Optional[str]) -> Optional[str]:
        """ The form LOCAL_REFERENCE_FIELD is stored in, so a lookup can be handed any other """
        return reference_id

    @property
    def short_identifier(self):
        """ What this record is known by - its local reference, else whatever external system
            manages it, else the pk. One identity rule for previews, search and node chips """
        local_reference = getattr(self, self.LOCAL_REFERENCE_FIELD, None)
        return local_reference or self.external_pk or f"({self.pk})"

    @property
    def external_manager(self):
        em = None
        if self.external_pk:
            em = self.external_pk.external_manager
        return em

    def can_write(self, _user) -> bool:
        cw = True
        if self.external_manager:
            cw &= self.external_manager.can_modify
        return cw


def patient_name(first_name, last_name):
    if not (first_name or last_name):
        return 'Unknown Patient'
    if first_name and last_name:
        return f"{first_name} {last_name}"
    return first_name or last_name


def patient_name_surname_first(first_name, last_name):
    if not (first_name or last_name):
        return 'Unknown Patient'
    if first_name and last_name:
        return f"{last_name}, {first_name}"
    return first_name or last_name

# PATIENT / SAMPLE AGES AND DATES:
# We recommend using dates rather than storing age/deceased as age will not keep
# up to date over time.
# Ie use Patient.date_of_birth and Patient.date_of_birth to store age/deceased
# and Specimen.collection_date - age and deceased can then be calculated from dates
#
# For client code, please use the properties age/deceased which hide the details
#
# But as this isn't always available - hardcode age etc using _underscore prefixed fields


# kwargs: patient - sent after Patient.merge() moves everything across
patient_merged_signal = Signal()


class Patient(GuardianPermissionsMixin, HasPhenotypeDescriptionMixin, ExternallyManagedModel, PreviewModelMixin):
    LOCAL_REFERENCE_FIELD = "patient_code"

    family_code = models.TextField(null=True, blank=True)
    patient_code = models.TextField(null=True, blank=True)
    first_name = models.TextField(null=True, blank=True)
    last_name = models.TextField(null=True, blank=True)
    date_of_birth = models.DateField(null=True, blank=True)
    date_of_death = models.DateField(null=True, blank=True)
    sex = models.CharField(max_length=1, choices=Sex.choices, default=Sex.UNKNOWN)
    phenotype = models.TextField(null=True, blank=True)
    affected = models.BooleanField(null=True)
    consanguineous = models.BooleanField(null=True)
    research_consent = models.BooleanField(null=True, blank=True)
    medicare = models.TextField(null=True, blank=True)
    billing_details = models.TextField(null=True, blank=True)

    # Address fields
    street_address = models.TextField(null=True, blank=True)
    suburb = models.TextField(null=True, blank=True)
    postcode = models.TextField(null=True, blank=True)  # Text for international - UK/Canada use alphanumeric
    state = models.TextField(null=True, blank=True)
    telephone = models.TextField(null=True, blank=True)  # Text to keep leading zeros and '+' country codes

    fake_data = models.ForeignKey(FakeData, null=True, on_delete=CASCADE)
    _deceased = models.BooleanField(null=True, blank=True)

    class Meta:
        # The least that identifies a patient (#2518) - Patient.match holds the import's stricter rule.
        # A blank or whitespace-only value counts as missing; isnull is needed as a NULL check passes
        constraints = [
            models.CheckConstraint(
                condition=Q(patient_code__isnull=False, patient_code__regex=r"\S")
                | Q(last_name__isnull=False, last_name__regex=r"\S"),
                name="patient_has_code_or_last_name",
                violation_error_message="A patient needs a patient code or last name",
            ),
        ]

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-user-injured"

    @classmethod
    def preview_if_url_visible(cls) -> Optional[str]:
        return 'patients'

    @property
    def preview(self) -> PreviewData:
        parts = []
        if self.sex:
            parts.append(PreviewKeyValue(key="Sex", value=self.sex.title()))
        if self.date_of_birth:
            parts.append(PreviewKeyValue(key="DOB", value=self.date_of_birth))
        if self.deceased:
            parts.append(PreviewKeyValue(value="deceased"))

        return self.preview_with(
            identifier=self.short_identifier,
            title=self.display_identity,
            summary_extra=parts
        )

    def can_write(self, user) -> bool:
        return ExternallyManagedModel.can_write(self, user) and GuardianPermissionsMixin.can_write(self, user)

    @classmethod
    def filter_writable_for_user(cls, user, pks=None):
        """ Batch can_write - a record an external manager owns is read only here """
        qs = super().filter_writable_for_user(user, pks=pks)
        return qs.exclude(external_pk__external_manager__can_modify=False)

    @classmethod
    def allow_group_permission_delete(cls) -> bool:
        # Deletable via the group_permissions delete view (can_write() still blocks externally-managed ones)
        return True

    def _merge_values(self) -> dict:
        """ attname -> value for each field that has been set (not empty, not its default) """
        values = {}
        for field in Patient._meta.concrete_fields:
            if field.primary_key or field.name in ("created", "modified"):
                continue
            value = getattr(self, field.attname)
            if value in (None, "") or (field.has_default() and value == field.get_default()):
                continue
            values[field.attname] = value
        return values

    def get_merge_conflicts(self, other: 'Patient') -> list[str]:
        """ Why other can't be merged into this patient - empty if it can """
        conflicts = []
        values = self._merge_values()
        other_values = other._merge_values()
        for attname in sorted(values.keys() & other_values.keys()):
            if values[attname] != other_values[attname]:
                conflicts.append(f"{attname} {values[attname]!r} vs {other_values[attname]!r}")

        # Specimen reference_id is unique per patient
        other_specimens = other.specimen_set.values_list("reference_id", flat=True)
        for reference_id in self.specimen_set.filter(reference_id__in=other_specimens) \
                .values_list("reference_id", flat=True):
            conflicts.append(f"both have specimen '{reference_id}'")
        return conflicts

    @transaction.atomic
    def merge(self, other: 'Patient', description: str, origin: str, user: Optional[User] = None):
        """ Folds other into this patient: everything pointing at it moves here, any field set only on
            other is copied, and other is deleted. Raises ValueError if get_merge_conflicts finds any """
        if conflicts := self.get_merge_conflicts(other):
            raise ValueError(f"Can't merge {other}({other.pk}) into {self}({self.pk}): {', '.join(conflicts)}")

        for rel in Patient._meta.related_objects:
            # One to one (PhenotypeDescription) goes with other, and is redone on save if phenotype is copied
            if rel.one_to_many:
                rel.related_model._base_manager.filter(**{rel.field.name: other}).update(**{rel.field.name: self})

        values = self._merge_values()
        copied = {attname: value for attname, value in other._merge_values().items() if attname not in values}
        other_description = f"{other}({other.pk})"
        self._move_object_permissions(other)
        other.delete()  # before copying, as external_pk is unique
        if copied:
            for attname, value in copied.items():
                setattr(self, attname, value)
            self.save(update_fields=list(copied))
        PatientModification.objects.create(patient=self, user=user, origin=origin,
                                           description=f"Merged {other_description}: {description}")
        patient_merged_signal.send(sender=Patient, patient=self)

    def _move_object_permissions(self, other: 'Patient'):
        """ Whoever could see other keeps seeing its specimens and samples, which now hang off this patient """
        content_type = ContentType.objects.get_for_model(Patient)
        for permission_model, holder in ((UserObjectPermission, "user"), (GroupObjectPermission, "group")):
            other_permissions = permission_model.objects.filter(content_type=content_type, object_pk=str(other.pk))
            for object_permission in other_permissions:
                permission_model.objects.get_or_create(content_type=content_type, object_pk=str(self.pk),
                                                       permission=object_permission.permission,
                                                       **{holder: getattr(object_permission, holder)})
            other_permissions.delete()

    @property
    def code(self):
        return self.patient_code or self.external_pk or f"Patient:{self.pk}"

    @property
    def name(self):
        return patient_name(self.first_name, self.last_name)

    @property
    def name_last_name_first(self):
        return patient_name_surname_first(self.first_name, self.last_name)

    @property
    def deceased(self):
        return any([self._deceased, self.date_of_death])

    @property
    def num_specimens(self):
        return self.specimen_set.count()

    @property
    def num_extractions(self):
        return Extraction.objects.filter(specimen__patient=self).count()

    @property
    def age(self):
        """ You may actually want to use Specimen.age_at_collection_date rather than this """
        return calculate_age(self.date_of_birth, self.date_of_death)

    @staticmethod
    def with_phenotype_text() -> QuerySet['Patient']:
        """ Patients whose phenotype text is matched - set, and free of settings.PATIENT_PHENOTYPE_EXCLUDE_STRING """
        patients = Patient.objects.filter(phenotype__isnull=False).exclude(phenotype='')
        if exclude_string := getattr(settings, "PATIENT_PHENOTYPE_EXCLUDE_STRING", None):
            patients = patients.exclude(phenotype__contains=exclude_string)
        return patients

    @staticmethod
    def for_ontology_term(user, ontology_term) -> QuerySet['Patient']:
        """ Patients the user can view whose phenotype text matched ontology_term """
        return Patient.filter_for_user(user).filter(**{PATIENT_ONTOLOGY_TERM_PATH: ontology_term}).order_by("id")

    def get_json_dict(self):
        args = {'sex': self.sex,
                'deceased': self.deceased}
        return args

    @staticmethod
    def match(first_name, last_name, sex=None, date_of_birth=None, user=None,
              patient_code=None) -> tuple['Patient', Optional[PatientRecordMatchType]]:
        """" A patient is identified by patient code alone, or by last name plus first name or date of birth (which
             must equal the patient's when there is no first name). A row with a code that no patient has only
             matches by name a patient without a code.
             We need to be able to match incomplete info, eg Male if provided matches Sex = M or Unknown.
             A match where the patient lacks a field the row has is PARTIAL """

        if patient_code:
            identity_q = Q(patient_code=patient_code)
        elif not last_name:
            msg = "Patient code or last name must be non-null!"
            raise ValueError(msg)
        elif not (first_name or date_of_birth):
            msg = "First name or date of birth is needed with the last name"
            raise ValueError(msg)
        else:
            identity_q = None

        if user:
            patients_queryset = Patient.filter_for_user(user)
        else:
            patients_queryset = Patient.objects.all()

        if identity_q is not None:
            patient = Patient._get_single_match(patients_queryset, identity_q)
            if patient:
                Patient._check_code_match_agrees(patient, first_name, last_name, date_of_birth)
                return patient, Patient._match_type(patient, sex, date_of_birth)
            if not last_name or not (first_name or date_of_birth):
                return None, None
            patients_queryset = patients_queryset.filter(Q(patient_code__isnull=True) | Q(patient_code=""))

        q = Q(last_name__iexact=last_name)
        if first_name:
            q &= Q(first_name__iexact=first_name)
            if date_of_birth is not None:
                q &= (Q(date_of_birth=date_of_birth) | Q(date_of_birth__isnull=True))
        else:
            q &= Q(date_of_birth=date_of_birth)
        if sex in Sex.FILLED_IN_CHOICES:
            q &= Sex.get_q_handling_unknown('sex', sex)

        patient = Patient._get_single_match(patients_queryset, q)
        if patient is None:
            return None, None
        return patient, Patient._match_type(patient, sex, date_of_birth)

    @staticmethod
    def _get_single_match(patients_queryset, q) -> Optional['Patient']:
        try:
            return patients_queryset.get(q)
        except Patient.DoesNotExist:
            return None
        except Patient.MultipleObjectsReturned as e:
            # Name them, so whoever is importing can fix the row or merge the patients
            matched = ", ".join(f"Patient:{p.pk} {p}" for p in patients_queryset.filter(q).order_by("pk"))
            raise Patient.MultipleObjectsReturned(f"Matched multiple patients: {matched}. "
                                                  "Add date of birth or sex to the row to tell them apart") from e

    @staticmethod
    def _check_code_match_agrees(patient, first_name, last_name, date_of_birth):
        """ A code matched, but a name or date of birth on both the row and the patient that differs means the
            code is on the wrong row (or the wrong patient) - refuse rather than attach it to someone else """
        differs = []
        for label, row_value, patient_value in [("first name", first_name, patient.first_name),
                                                ("last name", last_name, patient.last_name)]:
            if row_value and patient_value and row_value.upper() != patient_value.upper():
                differs.append(label)
        if isinstance(date_of_birth, str):
            date_of_birth = datetime.date.fromisoformat(date_of_birth)
        if date_of_birth and patient.date_of_birth and date_of_birth != patient.date_of_birth:
            differs.append("date of birth")
        if differs:
            msg = f"Patient code {patient.patient_code} matched Patient:{patient.pk} but the {', '.join(differs)} differs"
            raise ValueError(msg)

    @staticmethod
    def _match_type(patient, sex, date_of_birth) -> PatientRecordMatchType:
        if (date_of_birth and patient.date_of_birth is None) or (sex in Sex.FILLED_IN_CHOICES and patient.sex == Sex.UNKNOWN):
            return PatientRecordMatchType.PARTIAL
        return PatientRecordMatchType.EXACT

    def save(self, *args, **kwargs):
        ensure_mutally_exclusive_fields_not_set(self, "date_of_death", "_deceased")
        pheno_kwargs = HasPhenotypeDescriptionMixin.pop_kwargs(kwargs)
        super().save(*args, **kwargs)
        HasPhenotypeDescriptionMixin.save_phenotype(self, pheno_kwargs)

    def get_samples(self):
        """ Every sample that reaches this patient, either way round - the patient CSV sets patient
            and may leave extraction null, and a hand-set patient can disagree with the extraction's.
            Same union as SOURCE_LEVELS[PATIENT] in patients.sample_grouping;
            Sample is taken off the relation because snpdb imports this module. """
        sample_model = self.sample_set.model
        reaches_patient = Q(patient=self) | Q(extraction__specimen__patient=self)
        return sample_model.objects.filter(reaches_patient).distinct() \
            .select_related("vcf__genome_build", "extraction__specimen").order_by("vcf__date")

    @property
    def display_identity(self) -> str:
        """ The code alone when there is one - showing the name beside a de-identified code would
            re-identify the patient. Only a patient with no code is known by name """
        if code := self.patient_code or self.external_pk:
            return str(code)
        if self.first_name or self.last_name:
            return self.name_last_name_first
        return str(self.code)

    @staticmethod
    def display_identity_expression(prefix: str = "") -> Case:
        """ display_identity as SQL, so a grid can sort, filter and export on the identity its cells show.
            Change it alongside the property above - the two rules have to give the same answer.
            prefix: the lookup path to the patient, eg "patient__" from a model that has one """
        def field(name: str) -> str:
            return f"{prefix}{name}"

        def present(name: str) -> Q:
            """ Python treats a blank TextField as absent, so a NOT NULL test isn't enough """
            return Q(**{f"{field(name)}__isnull": False}) & ~Q(**{field(name): ""})

        has_first = present("first_name")
        has_last = present("last_name")
        return Case(
            # A prefixed path can arrive at no patient at all - that is a blank cell, not "Patient:"
            When(Q(**{f"{field('pk')}__isnull": True}), then=Value(None)),
            When(present("patient_code"), then=F(field("patient_code"))),
            When(Q(**{f"{field('external_pk')}__isnull": False}),
                 then=Concat(field("external_pk__code"), Value(" ("), field("external_pk__external_type"),
                             Value(")"), output_field=TextField())),
            When(has_first & has_last, then=Concat(field("last_name"), Value(", "), field("first_name"),
                                                   output_field=TextField())),
            When(has_first, then=F(field("first_name"))),
            When(has_last, then=F(field("last_name"))),
            default=Concat(Value("Patient:"), field("pk"), output_field=TextField()),
            output_field=TextField())

    def __str__(self):
        description = self.display_identity
        if self.sex != Sex.UNKNOWN:
            description += f" ({self.sex})"
        return description

    def get_absolute_url(self):
        return reverse('view_patient', kwargs={"patient_id": self.pk})

    @classmethod
    def get_listing_url(cls):
        return reverse('patients')


class PatientPopulation(TimeStampedModel):
    """ Can have many-to-one - e.g. Obama would have an entry for
        both AFRICAN_AFRICAN_AMERICAN and NON_FINNISH_EUROPEAN """
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    population = models.CharField(max_length=3, choices=PopulationGroup.choices)


class Tissue(TimeStampedModel):
    name = models.TextField()
    description = models.TextField()

    def __str__(self):
        return self.name


class NormalizedReferenceQuerySet(QuerySet):
    """ The bulk writes that skip save() still store reference_id in the model's canonical form """

    def bulk_create(self, objs, *args, **kwargs):
        for obj in objs:
            obj.reference_id = self.model.normalize_reference_id(obj.reference_id)
        return super().bulk_create(objs, *args, **kwargs)

    def update(self, **kwargs):
        if "reference_id" in kwargs:
            kwargs["reference_id"] = self.model.normalize_reference_id(kwargs["reference_id"])
        return super().update(**kwargs)


class AmbiguousExtraction(ValueError):
    """ A patient CSV row matches several of a specimen's extractions and has nothing to pick one with """

    def __init__(self, specimen, nucleic_acid_source, extractions):
        names = ", ".join(e.reference_id or f"unnamed:{e.pk}" for e in extractions)
        acid = f"{NucleicAcid(nucleic_acid_source).label} " if nucleic_acid_source else ""
        super().__init__(f"Specimen {specimen.reference_id} has {len(extractions)} {acid}extractions ({names}): "
                         "name one with Extraction Reference id or Extraction date")


class Specimen(GuardianPermissionsMixin, ExternallyManagedModel, PreviewModelMixin):
    """ Biological material collected from a patient - one tissue at one timepoint (block, blood draw).
        The nucleic acid taken off it lives on Extraction below """
    reference_id = models.TextField()
    description = models.TextField(null=True, blank=True)
    collected_by = models.TextField(null=True, blank=True)
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    tissue = models.ForeignKey(Tissue, null=True, blank=True, on_delete=SET_NULL)
    collection_date = models.DateTimeField(null=True, blank=True)
    received_date = models.DateTimeField(null=True, blank=True)
    # Per-specimen rather than per-tissue - the same tissue plays different roles in different tests
    tissue_status = models.CharField(max_length=1, choices=TissueStatus.choices, default=TissueStatus.UNKNOWN)
    # See note on patient / sample ages and dates above Patient model
    _age_at_collection_date = models.IntegerField(null=True, blank=True)

    class Meta:
        unique_together = ("patient", "reference_id")

    objects = NormalizedReferenceQuerySet.as_manager()

    @classmethod
    def normalize_reference_id(cls, reference_id: Optional[str]) -> Optional[str]:
        return normalize_specimen_reference(reference_id)

    def save(self, *args, **kwargs):
        self.reference_id = self.normalize_reference_id(self.reference_id)
        super().save(*args, **kwargs)

    @classmethod
    def get_permission_class(cls):
        return Patient

    def get_permission_object(self):
        # A specimen's confidentiality is its patient's
        return self.patient

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(patient__in=queryset)

    def can_write(self, user) -> bool:
        return ExternallyManagedModel.can_write(self, user) and GuardianPermissionsMixin.can_write(self, user)

    @classmethod
    def preview_icon(cls) -> str:
        # A blob against the extraction's stick - biological material vs the lab glassware taken off it
        return "fa-solid fa-droplet"

    @classmethod
    def preview_if_url_visible(cls) -> Optional[str]:
        return 'patients'

    @property
    def preview(self) -> PreviewData:
        parts = [PreviewKeyValue(key="Tissue status", value=self.get_tissue_status_display())]
        if self.collection_date:
            parts.append(PreviewKeyValue(key="Collected", value=self.collection_date))

        return self.preview_with(
            identifier=self.short_identifier,
            title=str(self.patient),
            summary_extra=parts
        )

    def get_absolute_url(self):
        return reverse('view_specimen', kwargs={"specimen_id": self.pk})

    @property
    def age_at_collection_date(self):
        if self._age_at_collection_date is not None:
            age = self._age_at_collection_date
        else:
            age = None
            # One of the dates should work....
            sample_date = self.collection_date or self.received_date
            if self.patient.date_of_birth and sample_date:
                age = calculate_age(self.patient.date_of_birth, sample_date)

        return age

    def get_or_create_extraction(self, nucleic_acid_source=None, sample=None, reference_id=None,
                                 extraction_date=None) -> 'Extraction':
        """ The extraction a patient CSV row describes, first hit wins:

            1. reference_id: the extraction with that reference on this specimen (the key the API upserts
               on), its acid and date set from the row where given - else a new one with all three.
               A reference beats the sample, so a row naming another arm relinks the sample
            2. The sample's current extraction, when it is on this specimen (acid updated in place) - so an
               existing file re-imports as before
            3. The specimen's extractions, filtered by the row's acid and date where given: one is it,
               several raise AmbiguousExtraction rather than guess (two DNA re-extractions look alike without
               a date). None, and the row gave an acid: an extraction with no acid yet is the row's to fill
               in. Otherwise a new one """
        extractions = self.extraction_set.order_by("pk")
        if reference_id := Extraction.normalize_reference_id(reference_id):
            if extraction := extractions.filter(reference_id=reference_id).first():
                # A blank cell leaves the stored value alone
                row_values = {"nucleic_acid_source": nucleic_acid_source, "extraction_date": extraction_date}
                changed = False
                for field, value in row_values.items():
                    if value is not None and getattr(extraction, field) != value:
                        setattr(extraction, field, value)
                        changed = True
                if changed:
                    extraction.save()
                return extraction
            return Extraction.objects.create(specimen=self, reference_id=reference_id,
                                             nucleic_acid_source=nucleic_acid_source,
                                             extraction_date=extraction_date)

        if sample and sample.extraction_id:
            if extraction := extractions.filter(pk=sample.extraction_id).first():
                if nucleic_acid_source and extraction.nucleic_acid_source != nucleic_acid_source:
                    extraction.nucleic_acid_source = nucleic_acid_source
                    extraction.save()
                return extraction

        candidates = extractions
        if nucleic_acid_source:
            candidates = candidates.filter(nucleic_acid_source=nucleic_acid_source)
        if extraction_date:
            if timezone.is_aware(extraction_date):
                extraction_date = timezone.localtime(extraction_date)
            candidates = candidates.filter(extraction_date__date=extraction_date.date())
        matches = list(candidates)
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise AmbiguousExtraction(self, nucleic_acid_source, matches)

        if nucleic_acid_source:
            # An extraction nobody has named an acid for is the row's to fill in; one already naming
            # a different acid is somebody else's arm
            unnamed_source = Q(nucleic_acid_source__isnull=True) | Q(nucleic_acid_source="")
            if extraction := extractions.filter(unnamed_source).first():
                extraction.nucleic_acid_source = nucleic_acid_source
                if extraction_date:
                    extraction.extraction_date = extraction_date
                extraction.save()
                return extraction

        return Extraction.objects.create(specimen=self, nucleic_acid_source=nucleic_acid_source,
                                         extraction_date=extraction_date)

    def __str__(self):
        s = self.reference_id
        if self.description:
            s += f" ({self.description})"
        return s


class Extraction(GuardianPermissionsMixin, ExternallyManagedModel, PreviewModelMixin):
    """ Nucleic acid taken off a Specimen - eg the DNA arm and the RNA arm of one tumour block.
        One extraction can be sequenced more than once (repeats, top-ups) so SequencingSample
        points here rather than the other way around """
    specimen = models.ForeignKey(Specimen, on_delete=CASCADE)
    # Local reference beside external_pk - not every deployment has a system to be managed by,
    # the same reason Patient carries patient_code and Specimen carries its own reference_id
    reference_id = models.TextField(null=True, blank=True)
    nucleic_acid_source = models.CharField(max_length=1, choices=NucleicAcid.choices, default=NucleicAcid.DNA, null=True, blank=True)
    extraction_date = models.DateTimeField(null=True, blank=True)

    class Meta:
        # Unnamed extractions are all distinct under Postgres, which is what the re-extraction
        # case wants - extraction_date tells those apart
        unique_together = ("specimen", "reference_id")

    objects = NormalizedReferenceQuerySet.as_manager()

    @classmethod
    def normalize_reference_id(cls, reference_id: Optional[str]) -> Optional[str]:
        return normalize_specimen_reference(reference_id)

    def save(self, *args, **kwargs):
        self.reference_id = self.normalize_reference_id(self.reference_id)
        super().save(*args, **kwargs)

    @classmethod
    def get_permission_class(cls):
        return Patient

    def get_permission_object(self):
        return self.specimen.patient

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(specimen__patient__in=queryset)

    def can_write(self, user) -> bool:
        return ExternallyManagedModel.can_write(self, user) and GuardianPermissionsMixin.can_write(self, user)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-flask"

    @classmethod
    def preview_if_url_visible(cls) -> Optional[str]:
        return 'patients'

    @property
    def preview(self) -> PreviewData:
        parts = []
        if self.nucleic_acid_source:
            parts.append(PreviewKeyValue(key="Nucleic acid", value=self.get_nucleic_acid_source_display()))
        if self.extraction_date:
            parts.append(PreviewKeyValue(key="Extracted", value=self.extraction_date))

        return self.preview_with(
            identifier=self.short_identifier,
            title=str(self.specimen.patient),
            summary_extra=parts
        )

    def get_absolute_url(self):
        return reverse('view_extraction', kwargs={"extraction_id": self.pk})

    def __str__(self):
        s = self.reference_id or str(self.specimen)
        if self.nucleic_acid_source:
            s += f" ({self.get_nucleic_acid_source_display()})"
        if self.extraction_date:
            s += f" {self.extraction_date.strftime(settings.DATE_FORMAT)}"
        return s


class ExtractionMatchMixin(models.Model):
    """ A claim about which Extraction a row belongs to, which may not be resolvable yet.

        Neither Sample nor SequencingSample is a TimeStampedModel, so extraction_match_date carries
        when the claim was parked - which is what promotes a stale Pending to Needs attention. """
    # Optional on a SequencingSample - seqauto isn't run everywhere, so Sample.extraction is the join
    # key. TODO: A sample may have >1 extractions (eg tumor/normal subtraction)
    extraction = models.ForeignKey(Extraction, null=True, blank=True, on_delete=SET_NULL)
    extraction_reference = models.JSONField(null=True, blank=True)
    extraction_match_status = models.CharField(max_length=1, choices=MatchStatus.choices,
                                               null=True, blank=True)
    extraction_match_error = models.TextField(null=True, blank=True)
    extraction_match_date = models.DateTimeField(null=True, blank=True)

    class Meta:
        abstract = True

    def apply_extraction_match(self, resolved: ResolvedReference, save=True) -> bool:
        """ Leaves a confirmed link alone, so a settled row never flaps back """
        if self.extraction_id:
            return False
        reference_json = resolved.reference.as_json()
        if self.extraction_match_date is None or reference_json != self.extraction_reference:
            # A new claim starts the clock - re-resolving the same one leaves it where it was, so an
            # unresolvable reference actually ages past the pending window rather than being renewed
            self.extraction_match_date = timezone.now()
        self.extraction_reference = reference_json
        self.extraction_match_status = resolved.status
        self.extraction_match_error = resolved.error
        if resolved.matched:
            self.extraction = resolved.obj
        if save:
            self.save()
        return True


def measure_value_description(value: Optional[float], unit: Optional[str], call: Optional[str]) -> str:
    """ A number and the lab's call, worded the one way - the specimen grid lists it after the
        measure's name, the case report build form shows it beside that measure's checkbox """
    description = ""
    if value is not None:
        unit = unit or ""
        separator = "" if unit in ("", "%") else " "
        description = f"{value}{separator}{unit}"
    if call:
        description = f"{description} ({call})".strip()
    return description


class SpecimenMeasure(GuardianPermissionsMixin, TimeStampedModel):
    """ A scalar measured on the material itself rather than by sequencing it - the pathologist's tumour
        content. A sequencing analysis' numbers are that analysis' own record
        (seqauto.models.DragenTSO500CombinedVariantOutput), since a re-analysis changes them.

        Written by whatever holds the lab's pathology (SA Path's Mocha sync), with the source's own
        record in source_payload so 'which report' stays answerable. """
    specimen = models.ForeignKey(Specimen, on_delete=CASCADE)
    measure_type = models.CharField(max_length=1, choices=SpecimenMeasureType.choices)
    value = models.FloatField(null=True, blank=True)
    unit = models.TextField(null=True, blank=True)      # eg '%'
    call = models.TextField(null=True, blank=True)      # the pathologist's call - 'No tumour'
    method = models.TextField(blank=True)               # where it came from - 'Mocha'
    source_payload = models.JSONField(default=dict, blank=True)  # the source's record, as sent
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    class Meta:
        # One current value per measure - the report wants a single tumour content, so a resend
        # replaces rather than accumulating
        unique_together = ("specimen", "measure_type")

    @classmethod
    def get_permission_class(cls):
        return Patient

    def get_permission_object(self):
        return self.specimen.patient

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(specimen__patient__in=queryset)

    @property
    def value_description(self) -> str:
        return measure_value_description(self.value, self.unit, self.call)

    def __str__(self):
        return " ".join(filter(None, [self.get_measure_type_display(), self.value_description]))


class PatientAttachment(TimeStampedModel):
    UPLOAD_PATH = 'patient_attachments'

    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    file = models.FileField(upload_to=UPLOAD_PATH,
                            storage=PrivateUploadStorage())
    file_type = models.CharField(max_length=1, choices=AttachmentFileType.CHOICES)
    thumbnail_path = models.TextField(null=True)  # This is created in patient_attachment_post_save

    def get_file_dict(self):
        basename = os.path.basename(self.file.path)
        image_url = reverse('view_patient_file_attachment', kwargs={'patient_attachment_id': self.pk})

        if self.file_type == AttachmentFileType.IMAGE:
            thumb_url = reverse('view_patient_file_attachment_thumbnail', kwargs={'patient_attachment_id': self.pk})
        else:
            thumb_url = None

        if os.path.exists(self.file.path):
            size = self.file.size
        else:
            size = 0

        return {
            'pk': self.pk,
            'name': basename,
            'size': size,
            'url': image_url,
            'thumbnailUrl': thumb_url,  # rendered as a poster thumbnail
            'deleteUrl': reverse('patient_file_delete', kwargs={'pk': self.pk}),
            'deleteType': 'POST',
        }


class PatientImport(TimeStampedModel):
    """ All the modifications hang off this. We don't want this to get deleted so modifications stay there
        Even after reload patient import """
    name = models.TextField()


class PatientModification(TimeStampedModel):
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    date = models.DateTimeField(auto_now_add=True)
    description = models.TextField(null=True)
    origin = models.CharField(max_length=1, choices=PatientRecordOriginType.CHOICES)
    patient_import = models.ForeignKey(PatientImport, null=True, on_delete=CASCADE)

    def get_import_url(self):
        url = None
        if self.patient_import:
            try:
                url = self.patient_import.patientrecords.get_absolute_url()
            except PatientRecords.DoesNotExist:
                pass
        return url

class PatientComment(TimeStampedModel):
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    comment = models.TextField()
    patient_modification = models.ForeignKey(PatientModification, on_delete=CASCADE)


class PatientColumns:
    PATIENT_FAMILY_CODE = 'Family Code'
    PATIENT_CODE = 'Patient Code (de-identified)'
    PATIENT_FIRST_NAME = 'Patient First Name'
    PATIENT_LAST_NAME = 'Patient Last Name (required)'
    DATE_OF_BIRTH = 'Date of Birth (DD-MM-YYYY)'
    DATE_OF_DEATH = 'Date of Death (DD-MM-YYYY)'
    DECEASED = 'Deceased (Y/N) (mutually exclusive to date of death)'
    SEX = 'Sex (F/M/Unknown)'
    AFFECTED = 'Affected (Y/N)'
    CONSANGUINEOUS = 'Consanguineous (Y/N)'
    PATIENT_PHENOTYPE = 'Patient Phenotype'
    SAMPLE_ID = 'Sample ID'
    SAMPLE_NAME = 'Sample Name (match on name)'
    SPECIMEN_REFERENCE_ID = 'Specimen Reference id'
    SPECIMEN_DESCRIPTION = 'Specimen Description'
    SPECIMEN_COLLECTED_BY = 'Specimen Collected by'
    SPECIMEN_COLLECTION_DATE = 'Specimen Collection date'
    SPECIMEN_RECEIVED_DATE = 'Specimen Received date'
    SPECIMEN_TISSUE_STATUS = 'Specimen Tissue status (Reference/Affected/Unknown)'
    SPECIMEN_AGE_AT_COLLECTION_DATE = 'Age at collection date (mutually exclusive to date of birth)'
    EXTRACTION_REFERENCE_ID = 'Extraction Reference id'
    EXTRACTION_NUCLEIC_ACID_SOURCE = 'Extraction Nucleic acid source (DNA/RNA)'
    EXTRACTION_DATE = 'Extraction date'

    COLUMN_DETAILS = [
        (PATIENT_FAMILY_CODE, "String", "Family Code"),
        (PATIENT_CODE, "String", "De-identified code safe to display in grids"),
        (PATIENT_FIRST_NAME, "String", "First Name to match/create a patient"),
        (PATIENT_LAST_NAME, "String", "Last Name to match/create a patient"),
        (DATE_OF_BIRTH, "Date", "Date to match/create a patient"),
        (DATE_OF_DEATH, "Date", ""),
        (DECEASED, "Deceased Y/N", "Only use this column if date of death is unavailable"),
        (SEX, "Sex", ""),
        (AFFECTED, 'Affected', 'Affected with conditions/phenotype'),
        (CONSANGUINEOUS, 'Consanguineous', 'Patient has parents with same recent ancestor'),
        (PATIENT_PHENOTYPE, "String", "Phenotype for patient"),
        (SAMPLE_ID, "Int", "ID to match samples"),
        (SAMPLE_NAME, "String", "Name to match samples"),
        (SPECIMEN_REFERENCE_ID, "String", ""),
        (SPECIMEN_DESCRIPTION, "String", ""),
        (SPECIMEN_COLLECTED_BY, "String", ""),
        (SPECIMEN_COLLECTION_DATE, "Date", ""),
        (SPECIMEN_RECEIVED_DATE, "Date", ""),
        (SPECIMEN_TISSUE_STATUS, "Tissue Status", "Role the material plays in the test"),
        (SPECIMEN_AGE_AT_COLLECTION_DATE, "Int", "Only use this column if date of birth is unavailable"),
        (EXTRACTION_REFERENCE_ID, "String", "Names one extraction of the specimen, eg the DNA and RNA arms of a "
                                            "block - one row per extraction, repeating the specimen columns. "
                                            "Unique within the specimen. Blank: the row's sample, else the nucleic "
                                            "acid and extraction date, say which extraction it is"),
        (EXTRACTION_NUCLEIC_ACID_SOURCE, "Nucleic Acid Source", ""),
        (EXTRACTION_DATE, "Date", "Tells two extractions of the same acid apart"),
    ]

    COLUMNS = [c[0] for c in COLUMN_DETAILS]

    SAMPLE_QUERYSET_PATH = {
        PATIENT_FAMILY_CODE: "patient__family_code",
        PATIENT_CODE: "patient__patient_code",
        PATIENT_FIRST_NAME: "patient__first_name",
        PATIENT_LAST_NAME: "patient__last_name",
        DATE_OF_BIRTH: "patient__date_of_birth",
        DATE_OF_DEATH: "patient__date_of_death",
        DECEASED: "patient___deceased",
        SEX: "patient__sex",
        AFFECTED: "patient__affected",
        CONSANGUINEOUS: "patient__consanguineous",
        PATIENT_PHENOTYPE: "patient__phenotype",
        SAMPLE_ID: "pk",
        SAMPLE_NAME: "name",
        SPECIMEN_REFERENCE_ID: "extraction__specimen__reference_id",
        SPECIMEN_DESCRIPTION: "extraction__specimen__description",
        SPECIMEN_COLLECTED_BY: "extraction__specimen__collected_by",
        SPECIMEN_COLLECTION_DATE: "extraction__specimen__collection_date",
        SPECIMEN_RECEIVED_DATE: "extraction__specimen__received_date",
        SPECIMEN_TISSUE_STATUS: "extraction__specimen__tissue_status",
        SPECIMEN_AGE_AT_COLLECTION_DATE: "extraction__specimen___age_at_collection_date",
        EXTRACTION_REFERENCE_ID: "extraction__reference_id",
        EXTRACTION_NUCLEIC_ACID_SOURCE: "extraction__nucleic_acid_source",
        EXTRACTION_DATE: "extraction__extraction_date",
    }


def single_string_to_first_last_name_q(single_string):
    """ Q on first_name/last_name for a free-text name, eg "Dr Jo Smith" - an initial also matches the first name """
    name = nameparser.HumanName(single_string)
    if not (name.first or name.last):
        raise ValueError(f"No name found in '{single_string}'")

    queries = []
    first_name_queries = []
    first_name = name.first
    if first_name:
        first_name_queries.append(Q(first_name__iexact=first_name))

    last_name = name.last
    if last_name:
        queries.append(Q(last_name__iexact=last_name))
        # If we have the last name, also accept an initial
        if first_name:
            first_name_queries.append(Q(first_name__istartswith=first_name[0]))

    if first_name_queries:
        queries.append(reduce(operator.or_, first_name_queries))
    return reduce(operator.and_, queries)


class Clinician(TimeStampedModel):
    email = models.TextField(null=True, blank=True)
    title = models.CharField(max_length=1, choices=Title.CHOICES, null=True, blank=True)
    first_name = models.TextField(null=True, blank=True)
    last_name = models.TextField(null=True)
    specialty = models.TextField(null=True, blank=True)
    phone = models.TextField(null=True, blank=True)
    user = models.OneToOneField(User, null=True, on_delete=SET_NULL)  # Set user to this, and the user is a clinician

    @staticmethod
    def match(clinician_string):
        """ Attempt to load by pulling first_name, last_name out of single string. Throws exception if not found """
        q = single_string_to_first_last_name_q(clinician_string)
        #logging.info(q)
        return Clinician.objects.get(q)

    @staticmethod
    def user_is_clinician(user):
        return user.is_authenticated and Clinician.objects.filter(user=user).exists()

    @property
    def name(self):
        title = None
        if self.title is not None:
            title = self.get_title_display()
        name_bits = [title, self.first_name, self.last_name]
        return ' '.join([x for x in name_bits if x])

    def __str__(self):
        return self.name

########
# CSV imports


class PatientRecords(TimeStampedModel):
    patient_import = models.OneToOneField(PatientImport, on_delete=CASCADE)

    def can_view(self, user) -> bool:
        return user.is_superuser or (self.user is not None and user == self.user)

    def check_can_view(self, user):
        if not self.can_view(user):
            msg = f"You do not have permission to access PatientRecords pk={self.pk}"
            raise PermissionDenied(msg)

    @property
    def user(self) -> Optional[User]:
        if file_upload := self.file_upload:
            return file_upload.user
        return None

    @property
    def file_upload(self):
        """ Records from old imports may have had their FileUpload removed """
        try:
            return self.uploadedpatientrecords.file_upload
        except ObjectDoesNotExist:
            return None

    def get_absolute_url(self):
        return reverse('view_patient_import', kwargs={"patient_records_id": self.pk})


class PatientRecord(TimeStampedModel):
    """ Record created from Uploaded Patient Records
        will be processed into Sample/Patients/Specimens etc. """
    patient_records = models.ForeignKey(PatientRecords, on_delete=CASCADE)
    record_id = models.IntegerField()  # From file
    valid = models.BooleanField(default=False)
    validation_message = models.TextField(blank=True, null=True)

    # Data from match/create
    sample = models.ForeignKey('snpdb.Sample', null=True, on_delete=SET_NULL)
    patient = models.ForeignKey(Patient, null=True, related_name='related_patient', on_delete=SET_NULL)
    patient_match = models.CharField(max_length=1, choices=PatientRecordMatchType.choices, null=True)
    specimen = models.ForeignKey(Specimen, null=True, related_name='related_specimen', on_delete=SET_NULL)
    specimen_match = models.CharField(max_length=1, choices=PatientRecordMatchType.choices, null=True)
    extraction = models.ForeignKey(Extraction, null=True, related_name='related_extraction', on_delete=SET_NULL)
    extraction_match = models.CharField(max_length=1, choices=PatientRecordMatchType.choices, null=True)

    # These are filled in from spreadsheet or manually entered
    sample_identifier = models.IntegerField(null=True)
    sample_name = models.TextField(null=True)
    patient_family_code = models.TextField(null=True)
    patient_code = models.TextField(null=True)
    patient_first_name = models.TextField(null=True)
    patient_last_name = models.TextField(null=True)
    date_of_birth = models.DateField(null=True)
    date_of_death = models.DateField(null=True)
    sex = models.CharField(max_length=1, choices=Sex.choices, null=True)
    affected = models.BooleanField(null=True, blank=True)
    consanguineous = models.BooleanField(null=True, blank=True)
    _deceased = models.BooleanField(null=True, blank=True)
    patient_phenotype = models.TextField(null=True)
    specimen_reference_id = models.TextField(null=True)
    specimen_description = models.TextField(null=True)
    specimen_collected_by = models.TextField(null=True)
    specimen_collection_date = models.TextField(null=True)
    specimen_received_date = models.TextField(null=True)
    specimen_tissue_status = models.CharField(max_length=1, choices=TissueStatus.choices, null=True)
    specimen_age_at_collection_date = models.IntegerField(null=True, blank=True)
    extraction_reference_id = models.TextField(null=True)
    extraction_nucleic_acid_source = models.CharField(max_length=1, choices=NucleicAcid.choices, null=True)
    extraction_date = models.TextField(null=True)  # As typed - Extraction.extraction_date holds the parsed value

    # Methods on this can include save / isvalid etc -
    # call out to helper funcs etc

    def get_absolute_url(self):
        return reverse('view_patient_record', kwargs={"pk": self.pk})


@receiver(models.signals.post_save, sender=PatientAttachment)
def patient_attachment_post_save(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    if created:
        instance.file_type = AttachmentFileType.get_type_for_file(instance.file.name)

        if instance.file_type == AttachmentFileType.IMAGE:
            options = {'size': (100, 100), 'crop': True}
            thumbnail_path = get_thumbnailer(instance.file).get_thumbnail(options).path
            instance.thumbnail_path = thumbnail_path
        instance.save()


class FollowLeadScientist(TimeStampedModel):
    follow = models.ForeignKey(User, related_name='followed', on_delete=CASCADE)
    user = models.ForeignKey(User, related_name='following', on_delete=CASCADE)


def get_lead_scientist_users_for_user(user):
    users_list = [user]
    ls_to_follow = FollowLeadScientist.objects.filter(user=user).values_list("follow", flat=True)
    users_list.extend(ls_to_follow)
    return users_list
