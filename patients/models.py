from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from django.db.models.deletion import CASCADE, SET_NULL
from django.dispatch.dispatcher import receiver
from django.urls.base import reverse
from django_extensions.db.models import TimeStampedModel
from easy_thumbnails.files import get_thumbnailer
import nameparser
import os

from annotation.models.has_phenotype_description_mixin import HasPhenotypeDescriptionMixin
from library.django_utils import single_string_to_first_last_name_q, \
    ensure_mutally_exclusive_fields_not_set
from library.django_utils.django_file_system_storage import PrivateUploadStorage
from library.django_utils.guardian_permissions_mixin import GuardianPermissionsMixin
from library.enums.file_attachments import AttachmentFileType
from library.enums.titles import Title
from library.utils import calculate_age
from patients.models_enums import NucleicAcid, Mutation, Sex, PopulationGroup

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


class ExternalPK(models.Model):
    code = models.TextField()
    external_type = models.TextField()
    external_manager = models.ForeignKey(ExternalModelManager, on_delete=CASCADE)

    class Meta:
        unique_together = ('code', 'external_type', 'external_manager')

    def __str__(self):
        return f"{self.code} ({self.external_type})"


class ExternallyManagedModel(TimeStampedModel):
    external_pk = models.OneToOneField(ExternalPK, null=True, on_delete=CASCADE)

    class Meta:
        abstract = True

    @property
    def external_manager(self):
        em = None
        if self.external_pk:
            em = self.external_pk.external_manager
        return em

    def can_write(self, user):
        cw = True
        if self.external_manager:
            cw &= self.external_manager.can_modify
        return cw


def patient_name(first_name, last_name):
    if not (first_name or last_name):
        return 'Unknown Patient'
    if first_name:
        return f"{first_name} {last_name}"
    return last_name

# PATIENT / SAMPLE AGES AND DATES:
# We recommend using dates rather than storing age/deceased as age will not keep
# up to date over time.
# Ie use Patient.date_of_birth and Patient.date_of_birth to store age/deceased
# and Specimen.collection_date - age and deceased can then be calculated from dates
#
# For client code, please use the properties age/deceased which hide the details
#
# But as this isn't always available - hardcode age etc using _underscore prefixed fields


class Patient(GuardianPermissionsMixin, HasPhenotypeDescriptionMixin, ExternallyManagedModel):
    family_code = models.TextField(null=True, blank=True)
    first_name = models.TextField(null=True, blank=True)
    last_name = models.TextField(null=True)
    date_of_birth = models.DateField(null=True, blank=True)
    date_of_death = models.DateField(null=True, blank=True)
    sex = models.CharField(max_length=1, choices=Sex.choices, default=Sex.UNKNOWN)
    phenotype = models.TextField(null=True, blank=True)
    affected = models.BooleanField(null=True)
    consanguineous = models.BooleanField(null=True)

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

    def can_write(self, user):
        return ExternallyManagedModel.can_write(self, user) and GuardianPermissionsMixin.can_write(self, user)

    @property
    def name(self):
        return patient_name(self.first_name, self.last_name)

    @property
    def deceased(self):
        return any([self._deceased, self.date_of_death])

    @property
    def condition_description(self):
        if self.date_of_death:
            return "dead (D.O.D. %s)" % self.date_of_death.strftime(settings.DATE_FORMAT)
        return "alive"

    @property
    def num_specimens(self):
        return self.specimen_set.count()

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
    def match(first_name, last_name, sex=None, date_of_birth=None, user=None):
        if not last_name:
            msg = "Last name must be non-null!"
            raise ValueError(msg)
        last_name = last_name.upper()
        if first_name:
            first_name = first_name.upper()
        patient_kwargs = {"last_name": last_name,
                          "first_name": first_name}
        if date_of_birth:
            patient_kwargs["date_of_birth"] = date_of_birth

        if sex in Sex.FILLED_IN_CHOICES:
            patient_kwargs['sex'] = sex

        try:
            if user:
                patients_queryset = Patient.filter_for_user(user)
            else:
                patients_queryset = Patient.objects.all()
            patient = patients_queryset.get(**patient_kwargs)
        except Patient.DoesNotExist:
            # TODO: Handle multiple patients - ie put in errors?
            patient = None

        return patient

    def save(self, **kwargs):
        ensure_mutally_exclusive_fields_not_set(self, "date_of_death", "_deceased")
        pheno_kwargs = HasPhenotypeDescriptionMixin.pop_kwargs(kwargs)
        super().save(**kwargs)
        HasPhenotypeDescriptionMixin.save_phenotype(self, pheno_kwargs)

    def get_samples(self):
        return self.sample_set.all().order_by("vcf__date")

    def __str__(self):
        description = self.name
        if self.sex != Sex.UNKNOWN:
            description += f" ({self.sex})"
        return description

    def get_absolute_url(self):
        return reverse('view_patient', kwargs={"patient_id": self.pk})

    @classmethod
    def get_listing_url(cls):
        return reverse('patients')


class PatientPopulation(models.Model):
    """ Can have many-to-one - eg Obama would have an entry for
        both AFRICAN_AFRICAN_AMERICAN and NON_FINNISH_EUROPEAN """
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    population = models.CharField(max_length=3, choices=PopulationGroup.choices)


class Tissue(models.Model):
    name = models.TextField()
    description = models.TextField()

    def __str__(self):
        return self.name


class Specimen(models.Model):
    """ Biological material in a test tube that was sequenced """
    reference_id = models.TextField(primary_key=True)
    description = models.TextField(null=True, blank=True)
    collected_by = models.TextField(null=True, blank=True)
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    tissue = models.ForeignKey(Tissue, null=True, blank=True, on_delete=SET_NULL)
    collection_date = models.DateTimeField(null=True, blank=True)
    received_date = models.DateTimeField(null=True, blank=True)
    mutation_type = models.CharField(max_length=1, choices=Mutation.choices, default=Mutation.GERMLINE, null=True, blank=True)
    nucleic_acid_source = models.CharField(max_length=1, choices=NucleicAcid.choices, default=NucleicAcid.DNA, null=True, blank=True)
    # See note on patient / sample ages and dates above Patient model
    _age_at_collection_date = models.IntegerField(null=True, blank=True)

    @property
    def age_at_collection_date(self):
        if self._age_at_collection_date:
            age = self._age_at_collection_date
        else:
            age = None
            # One of the dates should work....
            sample_date = self.collection_date or self.received_date
            if self.patient.date_of_birth and sample_date:
                age = calculate_age(self.patient.date_of_birth, sample_date)

        return age

    def __str__(self):
        return f"Specimen {self.reference_id}"


class PatientAttachment(models.Model):
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
            'name': basename,
            'size': size,
            'url': image_url,
            'thumbnailUrl': thumb_url,  # If thumbnail set, JFU displays in gallery
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


class PatientImport(models.Model):
    name = models.TextField()


class PatientModification(models.Model):
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    date = models.DateTimeField(auto_now_add=True)
    description = models.TextField(null=True)
    origin = models.CharField(max_length=1, choices=PatientRecordOriginType.CHOICES)
    patient_import = models.ForeignKey(PatientImport, null=True, on_delete=CASCADE)


class PatientComment(models.Model):
    patient = models.ForeignKey(Patient, on_delete=CASCADE)
    comment = models.TextField()
    patient_modification = models.ForeignKey(PatientModification, on_delete=CASCADE)


class PatientColumns:
    PATIENT_FAMILY_CODE = 'Family Code'
    PATIENT_FIRST_NAME = 'Patient First Name'
    PATIENT_LAST_NAME = 'Patient Last Name (or code) (required)'
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
    SPECIMEN_MUTATION_TYPE = 'Specimen Mutation type (Germline/Somatic)'
    SPECIMEN_NUCLEIC_ACID_SOURCE = 'Specimen Nucleic acid source (DNA/RNA)'
    SPECIMEN_AGE_AT_COLLECTION_DATE = 'Age at collection date (mutually exclusive to date of birth)'

    COLUMN_DETAILS = [
        (PATIENT_FAMILY_CODE, "String", "Family Code"),
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
        (SPECIMEN_MUTATION_TYPE, "Mutation Type", ""),
        (SPECIMEN_NUCLEIC_ACID_SOURCE, "Nucleic Acid Source", ""),
        (SPECIMEN_AGE_AT_COLLECTION_DATE, "Int", "Only use this column if date of birth is unavailable"),
    ]

    COLUMNS = [c[0] for c in COLUMN_DETAILS]

    SAMPLE_QUERYSET_PATH = {
        PATIENT_FAMILY_CODE: "patient__family_code",
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
        SPECIMEN_REFERENCE_ID: "specimen__reference_id",
        SPECIMEN_DESCRIPTION: "specimen__description",
        SPECIMEN_COLLECTED_BY: "specimen__collected_by",
        SPECIMEN_COLLECTION_DATE: "specimen__collection_date",
        SPECIMEN_RECEIVED_DATE: "specimen__received_date",
        SPECIMEN_MUTATION_TYPE: "specimen__mutation_type",
        SPECIMEN_NUCLEIC_ACID_SOURCE: "specimen__nucleic_acid_source",
        SPECIMEN_AGE_AT_COLLECTION_DATE: "specimen___age_at_collection_date",
    }


class Clinician(models.Model):
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
        except:
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


class PatientRecords(models.Model):
    patient_import = models.OneToOneField(PatientImport, on_delete=CASCADE)

    @property
    def user(self):
        return self.uploaded_file.user

    @property
    def uploaded_file(self):
        return self.uploadedpatientrecords.uploaded_file

    def get_absolute_url(self):
        return reverse('view_patient_records', kwargs={"patient_records_id": self.pk})


class PatientRecord(models.Model):
    """ Record created from Uploaded Patient Records
        will be processed into Sample/Patients/Specimens etc. """
    patient_records = models.ForeignKey(PatientRecords, on_delete=CASCADE)
    record_id = models.IntegerField()  # From file
    valid = models.BooleanField(default=False)
    validation_message = models.TextField(blank=True, null=True)

    # For each sample/patient/specimen there will be
    # EITHER a matched_ OR a created_ non null FK depending
    # On whether existing match or created with import
    matched_sample_id = models.IntegerField(null=True)  # Can't be a proper link due to cyclical dependencies
    matched_patient = models.ForeignKey(Patient, null=True, related_name='matched_patient', on_delete=CASCADE)
    matched_specimen = models.ForeignKey(Specimen, null=True, related_name='matched_specimen', on_delete=CASCADE)
    created_patient = models.ForeignKey(Patient, null=True, related_name='created_patient', on_delete=CASCADE)
    created_specimen = models.ForeignKey(Specimen, null=True, related_name='created_specimen', on_delete=CASCADE)

    # These are filled in from spreadsheet or manually entered
    sample_id = models.IntegerField(null=True)
    sample_name = models.TextField(null=True)
    patient_family_code = models.TextField(null=True)
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
    specimen_mutation_type = models.CharField(max_length=1, choices=Mutation.choices, null=True)
    specimen_nucleic_acid_source = models.CharField(max_length=1, choices=NucleicAcid.choices, null=True)
    specimen_age_at_collection_date = models.IntegerField(null=True, blank=True)

    # Methods on this can include save / isvalid etc -
    # call out to helper funcs etc


@receiver(models.signals.post_save, sender=PatientAttachment)
def patient_attachment_post_save(sender, instance, created, **kwargs):  # pylint: disable=unused-argument
    if created:
        instance.file_type = AttachmentFileType.get_type_for_file(instance.file.name)

        if instance.file_type == AttachmentFileType.IMAGE:
            options = {'size': (100, 100), 'crop': True}
            thumbnail_path = get_thumbnailer(instance.file).get_thumbnail(options).path
            instance.thumbnail_path = thumbnail_path
        instance.save()


class FollowLeadScientist(models.Model):
    follow = models.ForeignKey(User, related_name='followed', on_delete=CASCADE)
    user = models.ForeignKey(User, related_name='following', on_delete=CASCADE)


def get_lead_scientist_users_for_user(user):
    users_list = [user]
    ls_to_follow = FollowLeadScientist.objects.filter(user=user).values_list("follow")
    users_list.extend(ls_to_follow)
    return users_list
