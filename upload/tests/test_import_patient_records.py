"""upload/test_data/patient_upload.csv through file type detection and the import pipeline (#1745)"""
import os
from datetime import date

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, PatientRecord, Specimen
from patients.models_enums import NucleicAcid, Sex, TissueStatus
from snpdb.models import VCF, GenomeBuild, ImportSource, ImportStatus, ProcessingStatus, Sample
from upload.models import FileUpload, UploadedFileTypes, UploadPipeline
from upload.upload_processing import process_uploaded_file
from upload.uploaded_file_type import get_uploaded_file_type

PATIENT_UPLOAD_CSV = os.path.join(settings.BASE_DIR, "upload", "test_data", "patient_upload.csv")
SAMPLE_NAMES = ["FAM001_proband", "FAM001_mother", "FAM001_father", "ONC042_DNA", "ONC042_RNA"]


class TestPatientUploadCSV(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="patient_upload_user", password="x")
        cls.file_upload = FileUpload.objects.create(user=cls.user, name="patient_upload.csv", path=PATIENT_UPLOAD_CSV,
                                                    import_source=ImportSource.COMMAND_LINE)

    def test_recognised_as_patient_records(self):
        """ The csv extension is shared with gene lists, which claim any csv by default """
        self.assertEqual(UploadedFileTypes.PATIENT_RECORDS,
                         get_uploaded_file_type(self.file_upload, self.file_upload.name))


class TestPatientUploadImport(TestCase):
    maxDiff = None

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user(username="patient_import_user", password="x")
        vcf = VCF.objects.create(name="patient_upload.vcf", genotype_samples=len(SAMPLE_NAMES),
                                 genome_build=GenomeBuild.get_name_or_alias("GRCh37"),
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        assign_permission_to_user_and_groups(cls.user, vcf)
        for name in SAMPLE_NAMES:
            Sample.objects.create(name=name, vcf=vcf)
        file_upload = FileUpload.objects.create(user=cls.user, name="patient_upload.csv", path=PATIENT_UPLOAD_CSV,
                                                file_type=UploadedFileTypes.PATIENT_RECORDS,
                                                import_source=ImportSource.COMMAND_LINE)
        cls.upload_pipeline, _ = process_uploaded_file(file_upload, run_async=False)
        cls.upload_pipeline = UploadPipeline.objects.get(pk=cls.upload_pipeline.pk)

    def test_pipeline_succeeds(self):
        self.assertEqual(ProcessingStatus.SUCCESS, self.upload_pipeline.status, self.upload_pipeline.progress_status)
        self.assertEqual(7, self.upload_pipeline.items_processed)

    def test_every_record_valid(self):
        invalid = {r.record_id: r.validation_message for r in PatientRecord.objects.filter(valid=False)}
        self.assertEqual({}, invalid)

    def test_patients_created(self):
        self.assertEqual(6, Patient.objects.count())

    def test_dates_are_day_first(self):
        """ The header says DD-MM-YYYY, so 05/07/2015 is the 5th of July """
        junior = Patient.objects.get(first_name="JUNIOR", last_name="DOE")
        self.assertEqual(date(2015, 7, 5), junior.date_of_birth)

    def test_samples_matched_by_name(self):
        for sample in Sample.objects.filter(name__in=SAMPLE_NAMES):
            with self.subTest(sample=sample.name):
                self.assertIsNotNone(sample.patient)
                self.assertIsNotNone(sample.extraction)

    def test_one_specimen_two_arms(self):
        specimen = Specimen.objects.get(reference_id="SPEC-042")
        self.assertEqual(TissueStatus.AFFECTED, specimen.tissue_status)
        arms = {e.nucleic_acid_source for e in Extraction.objects.filter(specimen=specimen)}
        self.assertEqual({NucleicAcid.DNA, NucleicAcid.RNA}, arms)
        self.assertEqual(Sample.objects.get(name="ONC042_DNA").specimen, Sample.objects.get(name="ONC042_RNA").specimen)

    def test_tissue_status_from_header_values(self):
        self.assertEqual(TissueStatus.REFERENCE, Specimen.objects.get(reference_id="SPEC-002").tissue_status)
        self.assertEqual(TissueStatus.UNKNOWN, Specimen.objects.get(reference_id="SPEC-101").tissue_status)

    def test_deceased(self):
        pat = Patient.objects.get(last_name="PATIENTOVICH")
        self.assertEqual(date(2024, 1, 31), pat.date_of_death)
        self.assertTrue(pat.deceased)
        anon = Patient.objects.get(last_name="ANONYMOUS")
        self.assertIsNone(anon.date_of_death)
        self.assertTrue(anon.deceased)
        self.assertEqual(Sex.UNKNOWN, anon.sex)
        self.assertFalse(Patient.objects.get(first_name="JANE", last_name="DOE").deceased)

    def test_age_at_collection(self):
        self.assertEqual(45, Specimen.objects.get(reference_id="SPEC-101").age_at_collection_date)
