import os
from typing import Optional

import nameparser
from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.db import models
from django.db.models import Q
from django.db.models.deletion import CASCADE, SET_NULL
from django.dispatch.dispatcher import receiver
from django.urls.base import reverse
from django.utils import timezone
from django_extensions.db.models import TimeStampedModel

from annotation.models.has_phenotype_description_mixin import HasPhenotypeDescriptionMixin
from library.django_utils import (
    ensure_mutally_exclusive_fields_not_set,
    single_string_to_first_last_name_q,
)
from library.django_utils.django_file_system_storage import PrivateUploadStorage
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.enums.file_attachments import AttachmentFileType
from library.enums.titles import Title
from library.preview_request import PreviewData, PreviewKeyValue, PreviewModelMixin
from library.utils import calculate_age
from patients.external_references import ResolvedReference
from patients.models_enums import (
    MatchStatus,
    NucleicAcid,
    PatientRecordMatchType,
    PopulationGroup,
    Sex,
    SpecimenMeasureType,
    TissueStatus,
)

TEST_PATIENT_KWARGS = {"first_name": "PATIENT", "last_name": "TESTPATIENT"}


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
    telephone = models.IntegerField(null=True, blank=True)

    fake_data = models.ForeignKey(FakeData, null=True, on_delete=CASCADE)
    _deceased = models.BooleanField(null=True, blank=True)

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
            identifier=self.patient_code or self.external_pk or f"({self.pk})",
            title=self.name_last_name_first,
            summary_extra=parts
        )

    def can_write(self, user) -> bool:
        return ExternallyManagedModel.can_write(self, user) and GuardianPermissionsMixin.can_write(self, user)

    @classmethod
    def allow_group_permission_delete(cls) -> bool:
        # Deletable via the group_permissions delete view (can_write() still blocks externally-managed ones)
        return True

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
    def condition_description(self):
        if self.date_of_death:
            dod = self.date_of_death.strftime(settings.DATE_FORMAT)
            return f"dead (D.O.D. {dod})"
        if self._deceased:
            return "dead"
        return "alive"

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

    def _get_phenotype_input_text_field(self):
        # Implemented for HasPhenotypeDescriptionMixin
        return "phenotype"

    def _get_phenotype_description_relation_class_and_kwargs(self):
        # Implemented for HasPhenotypeDescriptionMixin
        # Stop circular import
        from annotation.models.models_phenotype_match import PatientTextPhenotype
        return PatientTextPhenotype, {"patient": self}

    def get_json_dict(self):
        args = {'sex': self.sex,
                'deceased': self.deceased}
        return args

    @staticmethod
    def match(first_name, last_name, sex=None, date_of_birth=None, user=None) -> tuple['Patient', Optional[PatientRecordMatchType]]:
        """" We need to be able to match incomplete info, eg Male if provided matches Sex = M or Unknown """

        if not last_name:
            msg = "Last name must be non-null!"
            raise ValueError(msg)
        q = Q(last_name__iexact=last_name)
        if first_name:
            q &= Q(first_name__iexact=first_name)
        if date_of_birth is not None:
            q &= (Q(date_of_birth=date_of_birth) | Q(date_of_birth__isnull=True))
        if sex in Sex.FILLED_IN_CHOICES:
            q &= Sex.get_q_handling_unknown('sex', sex)

        try:
            if user:
                patients_queryset = Patient.filter_for_user(user)
            else:
                patients_queryset = Patient.objects.all()

            patient = patients_queryset.get(q)
            if (date_of_birth and patient.date_of_birth is None) or (sex in Sex.FILLED_IN_CHOICES and patient.sex == Sex.UNKNOWN):
                match_type = PatientRecordMatchType.PARTIAL
            else:
                match_type = PatientRecordMatchType.EXACT
        except Patient.DoesNotExist:
            # TODO: Handle multiple patients - ie put in errors?
            patient = None
            match_type = None
        return patient, match_type

    def save(self, *args, **kwargs):
        ensure_mutally_exclusive_fields_not_set(self, "date_of_death", "_deceased")
        pheno_kwargs = HasPhenotypeDescriptionMixin.pop_kwargs(kwargs)
        super().save(*args, **kwargs)
        HasPhenotypeDescriptionMixin.save_phenotype(self, pheno_kwargs)

    def get_samples(self):
        return self.sample_set.all().select_related("vcf", "extraction__specimen").order_by("vcf__date")

    def __str__(self):
        # De-identified patients have no name, so fall back to the code they're known by
        if self.first_name or self.last_name:
            description = self.name
        else:
            description = str(self.code)
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
        return "fa-solid fa-vial"

    @classmethod
    def preview_if_url_visible(cls) -> Optional[str]:
        return 'patients'

    @property
    def preview(self) -> PreviewData:
        parts = [PreviewKeyValue(key="Tissue status", value=self.get_tissue_status_display())]
        if self.collection_date:
            parts.append(PreviewKeyValue(key="Collected", value=self.collection_date))

        return self.preview_with(
            identifier=self.reference_id or self.external_pk or f"({self.pk})",
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

    def get_or_create_extraction(self, nucleic_acid_source=None, sample=None) -> 'Extraction':
        """ The patient CSV describes a single extraction per specimen - the one the row's sample came
            off. Re-importing that row corrects that extraction in place, so the sample is what says
            which extraction the row is talking about.

            Without it, a row for a second sample would repurpose the first sample's extraction, and
            both arms of a TSO 500 specimen would end up on one. """
        extractions = self.extraction_set.order_by("pk")
        if sample and sample.extraction_id:
            if extraction := extractions.filter(pk=sample.extraction_id).first():
                if nucleic_acid_source and extraction.nucleic_acid_source != nucleic_acid_source:
                    extraction.nucleic_acid_source = nucleic_acid_source
                    extraction.save()
                return extraction

        if extraction := extractions.filter(nucleic_acid_source=nucleic_acid_source).first():
            return extraction

        if not nucleic_acid_source:
            # The CSV didn't say, so whatever is already there answers for it
            if extraction := extractions.first():
                return extraction
        else:
            # An extraction nobody has named an acid for is the CSV's to fill in; one already naming
            # a different acid is somebody else's arm
            unnamed_source = Q(nucleic_acid_source__isnull=True) | Q(nucleic_acid_source="")
            if extraction := extractions.filter(unnamed_source).first():
                extraction.nucleic_acid_source = nucleic_acid_source
                extraction.save()
                return extraction

        return Extraction.objects.create(specimen=self, nucleic_acid_source=nucleic_acid_source)

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
            identifier=self.reference_id or self.external_pk or f"({self.pk})",
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


class SpecimenMeasure(GuardianPermissionsMixin, TimeStampedModel):
    """ A scalar measured on the material rather than on any one variant - TMB, MSI, GIS.

        Both the score and the call are stored. HRD gives a genomic instability score and no
        positive/negative call: the threshold that turns it into one is lab policy, not vendor output,
        so without the raw value plus the threshold applied and by whom a later re-interpretation is
        unreconstructable. Mirrors somatic:tmb_value beside somatic:tmb_status in the evidence keys. """
    specimen = models.ForeignKey(Specimen, on_delete=CASCADE)
    # Which arm produced the number - enrichment, since the measure describes the specimen
    extraction = models.ForeignKey(Extraction, null=True, blank=True, on_delete=SET_NULL)
    measure_type = models.CharField(max_length=1, choices=SpecimenMeasureType.choices)
    value = models.FloatField(null=True, blank=True)
    unit = models.TextField(null=True, blank=True)      # eg 'mut/Mb', '%'
    call = models.TextField(null=True, blank=True)      # the lab's call - 'High', 'Stable'
    threshold = models.TextField(null=True, blank=True)  # the threshold that produced the call
    threshold_source = models.TextField(null=True, blank=True)  # whose policy set it
    method = models.TextField(blank=True)               # tool and version the client transcribed from
    source_payload = models.JSONField(default=dict, blank=True)  # raw, so 'which file' stays answerable
    measured_date = models.DateTimeField(null=True, blank=True)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)

    class Meta:
        # One current value per measure - the report pulls these together and wants a single MSI, so a
        # resend replaces rather than accumulating. Not keyed on the nullable extraction FK: Postgres
        # treats nulls as distinct, so that would let a resend duplicate instead of updating
        unique_together = ("specimen", "measure_type")

    @classmethod
    def get_permission_class(cls):
        return Patient

    def get_permission_object(self):
        return self.specimen.patient

    @classmethod
    def _filter_from_permission_object_qs(cls, queryset):
        return cls.objects.filter(specimen__patient__in=queryset)

    def __str__(self):
        description = self.get_measure_type_display()
        if self.value is not None:
            description += f" {self.value}{self.unit or ''}"
        if self.call:
            description += f" ({self.call})"
        return description


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


class PatientRecordOriginType:
    UPLOADED_CSV = 'C'
    MANUAL_VG_GUI = 'G'
    EXTERNAL_DATABASE = 'E'
    INTERNAL_VG = 'I'

    CHOICES = (
        (UPLOADED_CSV, 'Uploaded CSV'),
        (INTERNAL_VG, 'Internal change'),
        (MANUAL_VG_GUI, 'Manual change by user'),
        (EXTERNAL_DATABASE, 'External Database'),
    )


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
    SPECIMEN_NUCLEIC_ACID_SOURCE = 'Specimen Nucleic acid source (DNA/RNA)'
    SPECIMEN_AGE_AT_COLLECTION_DATE = 'Age at collection date (mutually exclusive to date of birth)'

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
        (SPECIMEN_NUCLEIC_ACID_SOURCE, "Nucleic Acid Source", ""),
        (SPECIMEN_AGE_AT_COLLECTION_DATE, "Int", "Only use this column if date of birth is unavailable"),
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
        SPECIMEN_NUCLEIC_ACID_SOURCE: "extraction__nucleic_acid_source",
        SPECIMEN_AGE_AT_COLLECTION_DATE: "extraction__specimen___age_at_collection_date",
    }


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
    def cleaned_get_or_create(clinician_string):
        try:
            clinician = Clinician.match(clinician_string)
        except Clinician.DoesNotExist:
            kwargs = {}
            name = nameparser.HumanName(clinician_string)
            if name.title:
                title = {v: k for k, v in Title.CHOICES}.get(name.title)
                if title:
                    kwargs["title"] = title

            first_name = name.first
            if name.middle:
                first_name += name.middle
            kwargs["first_name"] = first_name.upper()

            last_name = name.last
            if last_name:
                kwargs["last_name"] = last_name.upper()

            clinician = Clinician.objects.create(**kwargs)

        return clinician

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

    # These are filled in from spreadsheet or manually entered
    sample_identifier = models.IntegerField(null=True)
    sample_name = models.TextField(null=True)
    patient_family_code = models.TextField(null=True)
    patient_code = models.TextField(null=True)
    patient_first_name = models.TextField(null=True)
    patient_last_name = models.TextField()
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
    specimen_nucleic_acid_source = models.CharField(max_length=1, choices=NucleicAcid.choices, null=True)
    specimen_age_at_collection_date = models.IntegerField(null=True, blank=True)

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
            from easy_thumbnails.files import get_thumbnailer
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
