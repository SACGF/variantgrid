"""
Specimen -> Extraction -> Sample (#1704)
"""
import os

from django.contrib.auth.models import User
from django.test import TestCase
from django.utils import timezone

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.import_records import process_record
from patients.models import (
    Extraction,
    Patient,
    PatientColumns,
    PatientImport,
    PatientRecords,
    Specimen,
)
from patients.models_enums import Mutation, NucleicAcid
from patients.views_autocomplete import ExtractionAutocompleteView, SpecimenAutocompleteView
from snpdb.models import VCF, GenomeBuild, ImportSource, ImportStatus, Sample
from upload.models import FileUpload, UploadedFileTypes, UploadedPatientRecords

_FAKE_CSV = os.path.join(os.path.dirname(__file__), "test_data", "fake_patient_records.csv")


def _make_row(**overrides):
    row = dict.fromkeys(PatientColumns.COLUMNS)
    row[PatientColumns.PATIENT_LAST_NAME] = "EXTRACTLAST"
    row[PatientColumns.PATIENT_FIRST_NAME] = "EXTRACTFIRST"
    row.update(overrides)
    return row


class TestSpecimenExtraction(TestCase):
    """ One specimen, a DNA arm and an RNA arm - the TSO 500 shape that Specimen.nucleic_acid_source
        used to force into two unrelated specimens """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("extraction_user", password="x")
        cls.patient = Patient.objects.create(first_name="DNARNA", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000001", patient=cls.patient,
                                               mutation_type=Mutation.SOMATIC)
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        cls.rna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001B",
                                            nucleic_acid_source=NucleicAcid.RNA)

        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        vcf = VCF.objects.create(name="extraction_test.vcf", genotype_samples=2, genome_build=genome_build,
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        cls.dna_sample = Sample.objects.create(name="dna_arm", vcf=vcf, patient=cls.patient, extraction=cls.dna)
        cls.rna_sample = Sample.objects.create(name="rna_arm", vcf=vcf, patient=cls.patient, extraction=cls.rna)

    def test_both_arms_reach_one_specimen(self):
        self.assertEqual(self.dna_sample.extraction.specimen, self.specimen)
        self.assertEqual(self.rna_sample.extraction.specimen, self.specimen)
        self.assertEqual(self.dna_sample.specimen, self.rna_sample.specimen)

    def test_specimen_holds_dna_and_rna(self):
        sources = set(self.specimen.extraction_set.values_list("nucleic_acid_source", flat=True))
        self.assertEqual(sources, {NucleicAcid.DNA, NucleicAcid.RNA})

    def test_samples_selectable_by_specimen(self):
        sample_names = set(Sample.objects.filter(extraction__specimen=self.specimen).values_list("name", flat=True))
        self.assertEqual(sample_names, {"dna_arm", "rna_arm"})

    def test_reference_id_unique_per_patient_not_globally(self):
        """ R4 - the surrogate PK is what lets 2 patients use the same reference """
        other_patient = Patient.objects.create(first_name="OTHER", last_name="PATIENT")
        other_specimen = Specimen.objects.create(reference_id="2600000001", patient=other_patient)
        self.assertNotEqual(other_specimen.pk, self.specimen.pk)

    def test_get_or_create_extraction_matches_nucleic_acid_source(self):
        """ A specimen with both arms must not have one of them repurposed """
        self.assertEqual(self.specimen.get_or_create_extraction(NucleicAcid.RNA), self.rna)
        self.assertEqual(self.specimen.extraction_set.count(), 2)


class TestPatientCSVExtractionRoundTrip(TestCase):
    """ R6 - the existing SPECIMEN_NUCLEIC_ACID_SOURCE column implies a single extraction """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("csv_extraction_user", password="x")
        cls.patient = Patient.objects.create(first_name="EXTRACTFIRST", last_name="EXTRACTLAST")
        assign_permission_to_user_and_groups(cls.user, cls.patient)

        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        vcf = VCF.objects.create(name="csv_extraction_test.vcf", genotype_samples=1, genome_build=genome_build,
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        cls.sample = Sample.objects.create(name="csv_sample", vcf=vcf, import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(cls.user, cls.sample)

    def setUp(self):
        patient_import = PatientImport.objects.create(name="extraction_import")
        self.patient_records = PatientRecords.objects.create(patient_import=patient_import)
        file_upload = FileUpload.objects.create(user=self.user, name="extraction_import_file", path=_FAKE_CSV,
                                                file_type=UploadedFileTypes.PATIENT_RECORDS,
                                                import_source=ImportSource.COMMAND_LINE)
        UploadedPatientRecords.objects.create(file_upload=file_upload, patient_records=self.patient_records)

    def _import_row(self, **overrides):
        row = _make_row(**{
            PatientColumns.SAMPLE_NAME: "csv_sample",
            PatientColumns.SPECIMEN_REFERENCE_ID: "CSVSPEC001",
            PatientColumns.SPECIMEN_DESCRIPTION: "left femur",
            PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE: "DNA",
            **overrides,
        })
        process_record(self.patient_records, record_id=1, row=row)

    def test_import_creates_extraction_under_specimen(self):
        self._import_row()

        specimen = Specimen.objects.get(reference_id="CSVSPEC001")
        self.assertEqual(specimen.patient, self.patient)
        self.assertEqual(specimen.description, "left femur")

        extraction = Extraction.objects.get(specimen=specimen)
        self.assertEqual(extraction.nucleic_acid_source, NucleicAcid.DNA)

    def test_sample_reaches_specimen_via_extraction(self):
        self._import_row()

        self.sample.refresh_from_db()
        self.assertEqual(self.sample.extraction.specimen.reference_id, "CSVSPEC001")
        self.assertEqual(self.sample.specimen.reference_id, "CSVSPEC001")

    def test_reimport_reuses_the_same_extraction(self):
        self._import_row()
        extraction_pk = Extraction.objects.get(specimen__reference_id="CSVSPEC001").pk

        self._import_row()
        self.assertEqual(Extraction.objects.filter(specimen__reference_id="CSVSPEC001").count(), 1)
        self.assertEqual(Extraction.objects.get(specimen__reference_id="CSVSPEC001").pk, extraction_pk)

    def test_reimport_updates_nucleic_acid_source_in_place(self):
        self._import_row()
        self._import_row(**{PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE: "RNA"})

        extractions = Extraction.objects.filter(specimen__reference_id="CSVSPEC001")
        self.assertEqual(extractions.count(), 1)
        self.assertEqual(extractions.get().nucleic_acid_source, NucleicAcid.RNA)

    def test_blank_nucleic_acid_source_still_links_sample(self):
        self._import_row(**{PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE: None})

        self.sample.refresh_from_db()
        self.assertEqual(self.sample.extraction.specimen.reference_id, "CSVSPEC001")
        self.assertIsNone(self.sample.extraction.nucleic_acid_source)

    def test_export_columns_round_trip(self):
        """ The CSV written back out reads the same values through the new relations """
        self._import_row()

        values = Sample.objects.filter(pk=self.sample.pk).values(*PatientColumns.SAMPLE_QUERYSET_PATH.values()).get()
        self.assertEqual(values[PatientColumns.SAMPLE_QUERYSET_PATH[PatientColumns.SPECIMEN_REFERENCE_ID]],
                         "CSVSPEC001")
        self.assertEqual(values[PatientColumns.SAMPLE_QUERYSET_PATH[PatientColumns.SPECIMEN_DESCRIPTION]],
                         "left femur")
        self.assertEqual(values[PatientColumns.SAMPLE_QUERYSET_PATH[PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE]],
                         NucleicAcid.DNA)


class TestAutocompleteForwarding(TestCase):
    """ Picking one level of Patient -> Specimen -> Extraction narrows the others """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("autocomplete_user", password="x")
        cls.patient = Patient.objects.create(first_name="FORWARD", last_name="PATIENT")
        cls.other_patient = Patient.objects.create(first_name="OTHER", last_name="PATIENT")
        for patient in (cls.patient, cls.other_patient):
            assign_permission_to_user_and_groups(cls.user, patient)

        cls.specimen = Specimen.objects.create(reference_id="2600000001", patient=cls.patient)
        cls.other_specimen = Specimen.objects.create(reference_id="2600000002", patient=cls.patient)
        cls.other_patient_specimen = Specimen.objects.create(reference_id="2600000003",
                                                             patient=cls.other_patient)

        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        cls.rna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000001B",
                                            nucleic_acid_source=NucleicAcid.RNA)
        cls.other_extraction = Extraction.objects.create(specimen=cls.other_specimen,
                                                         reference_id="2600000002C")
        Extraction.objects.create(specimen=cls.other_patient_specimen, reference_id="2600000003C")

    def _queryset(self, view_class, **forwarded):
        view = view_class()
        view.forwarded = forwarded
        return view.get_user_queryset(self.user)

    def test_extraction_forwards_patient(self):
        extractions = self._queryset(ExtractionAutocompleteView, patient=self.patient.pk)
        self.assertEqual(set(extractions), {self.dna, self.rna, self.other_extraction})

    def test_extraction_forwards_specimen(self):
        """ Both arms of the forwarded specimen, and nothing from its sibling """
        extractions = self._queryset(ExtractionAutocompleteView, specimen=self.specimen.pk)
        self.assertEqual(set(extractions), {self.dna, self.rna})

    def test_specimen_forwards_extraction(self):
        specimens = self._queryset(SpecimenAutocompleteView, extraction=self.rna.pk)
        self.assertEqual(set(specimens), {self.specimen})

    def test_specimen_forwards_patient(self):
        specimens = self._queryset(SpecimenAutocompleteView, patient=self.patient.pk)
        self.assertEqual(set(specimens), {self.specimen, self.other_specimen})

    def test_nothing_forwarded_returns_all_visible(self):
        self.assertEqual(self._queryset(SpecimenAutocompleteView).count(), 3)
        self.assertEqual(self._queryset(ExtractionAutocompleteView).count(), 4)

    def test_extraction_searchable_by_its_own_reference(self):
        """ The TSO 500 container suffix lives on the extraction, not the specimen """
        self.assertIn('reference_id', ExtractionAutocompleteView.fields)
