"""upload/test_data/patient_upload.csv through file type detection and the import pipeline (#1745)"""
import os
import tempfile
from datetime import date
from unittest.mock import patch

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.enums import NucleicAcid, Sex, TissueStatus
from patients.import_records import get_patient_record_imports_dataframe
from patients.models import Extraction, Patient, PatientColumns, PatientRecord, Specimen
from snpdb.enums import ImportSource, ImportStatus, ProcessingStatus
from snpdb.models import VCF, GenomeBuild, Sample
from upload.enums import UploadedFileTypes
from upload.models import FileUpload, UploadPipeline
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

    def test_legacy_columns(self):
        """ A pre-#1870 file: nucleic acid under its Specimen name, no extraction reference or date """
        df = get_patient_record_imports_dataframe(PATIENT_UPLOAD_CSV)
        legacy_df = df.drop(columns=[PatientColumns.EXTRACTION_REFERENCE_ID, PatientColumns.EXTRACTION_DATE])
        legacy_df = legacy_df.rename(columns={PatientColumns.EXTRACTION_NUCLEIC_ACID_SOURCE:
                                              'Specimen Nucleic acid source (DNA/RNA)'})
        with tempfile.NamedTemporaryFile(suffix=".csv") as legacy_csv:
            legacy_df.to_csv(legacy_csv.name, index=False)
            imported_df = get_patient_record_imports_dataframe(legacy_csv.name)
        self.assertEqual(set(PatientColumns.COLUMNS), set(imported_df.columns))
        self.assertEqual(list(df[PatientColumns.EXTRACTION_NUCLEIC_ACID_SOURCE]),
                         list(imported_df[PatientColumns.EXTRACTION_NUCLEIC_ACID_SOURCE]))
        self.assertTrue(imported_df[PatientColumns.EXTRACTION_REFERENCE_ID].isna().all())

    def test_whitespace_cell_is_blank(self):
        """ A cell of only spaces reads as None, the same as an empty one """
        df = get_patient_record_imports_dataframe(PATIENT_UPLOAD_CSV)
        df[PatientColumns.PATIENT_LAST_NAME] = "  "
        with tempfile.NamedTemporaryFile(suffix=".csv") as spaces_csv:
            df.to_csv(spaces_csv.name, index=False)
            imported_df = get_patient_record_imports_dataframe(spaces_csv.name)
        self.assertTrue(all(v is None for v in imported_df[PatientColumns.PATIENT_LAST_NAME]))


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
        with patch("patients.import_records.reconcile_pending_extractions.delay") as reconcile:
            cls.upload_pipeline, _ = process_uploaded_file(file_upload, run_async=False)
        cls.reconcile_calls = reconcile.call_count
        cls.upload_pipeline = UploadPipeline.objects.get(pk=cls.upload_pipeline.pk)

    def test_pipeline_succeeds(self):
        self.assertEqual(ProcessingStatus.SUCCESS, self.upload_pipeline.status, self.upload_pipeline.progress_status)
        self.assertEqual(8, self.upload_pipeline.items_processed)

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
        """ Named arms plus a sample-less re-extraction, told apart by their Extraction Reference id """
        specimen = Specimen.objects.get(reference_id="2600000001")
        self.assertEqual(TissueStatus.AFFECTED, specimen.tissue_status)
        arms = {e.reference_id: (e.nucleic_acid_source, timezone.localdate(e.extraction_date))
                for e in Extraction.objects.filter(specimen=specimen)}
        self.assertEqual({"2600000001C": (NucleicAcid.DNA, date(2025, 5, 6)),
                          "2600000001B": (NucleicAcid.RNA, date(2025, 5, 6)),
                          "2600000001D": (NucleicAcid.DNA, date(2025, 6, 20))}, arms)
        self.assertEqual("2600000001C", Sample.objects.get(name="ONC042_DNA").extraction.reference_id)
        self.assertEqual("2600000001B", Sample.objects.get(name="ONC042_RNA").extraction.reference_id)

    def test_new_extractions_reconcile_parked_samples(self):
        self.assertEqual(1, self.reconcile_calls)

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
