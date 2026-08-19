"""
Specimen -> Extraction -> Sample (#1704), and Specimen.tissue_status (variantgrid_private#2447)
"""
import os
from datetime import datetime

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse
from django.utils import timezone

from classification.autopopulate_evidence_keys.evidence_from_sample_and_patient import (
    get_evidence_fields_for_sample_and_patient,
)
from classification.enums import SpecialEKeys
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.forms import ExtractionForm
from patients.import_records import process_record
from patients.models import (
    Extraction,
    ExternalModelManager,
    ExternalPK,
    Patient,
    PatientColumns,
    PatientImport,
    PatientRecord,
    PatientRecords,
    Specimen,
)
from patients.models_enums import NucleicAcid, TissueStatus
from patients.views import CREATE_EXTRACTION
from patients.views_autocomplete import ExtractionAutocompleteView, SpecimenAutocompleteView
from snpdb.models import VCF, GenomeBuild, ImportSource, ImportStatus, Sample
from snpdb.search import SearchResultMatchStrength, search_data
from snpdb.models.models_enums import VariantsType
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
                                               tissue_status=TissueStatus.AFFECTED)
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

    def test_extraction_date_distinguishes_a_re_extraction(self):
        """ Two DNA extractions off one specimen are told apart by when they were taken """
        redo = Extraction.objects.create(specimen=self.specimen, nucleic_acid_source=NucleicAcid.DNA,
                                         extraction_date=timezone.make_aware(datetime(2026, 3, 14)))
        self.dna.extraction_date = timezone.make_aware(datetime(2026, 1, 20))
        self.dna.save()

        self.assertNotEqual(str(redo), str(self.dna))
        dna_qs = self.specimen.extraction_set.filter(nucleic_acid_source=NucleicAcid.DNA)
        self.assertEqual(list(dna_qs.order_by("-extraction_date")), [redo, self.dna])

    def test_local_reference_and_external_pk_coexist(self):
        """ external_pk is nullable because not every deployment has a system managing it, so each
            level carries a local reference too """
        external_manager = ExternalModelManager.objects.create(name="test_lims")
        external_pk = ExternalPK.objects.create(code="LIMS-EXT-1", external_type="extraction",
                                                external_manager=external_manager)
        self.dna.external_pk = external_pk
        self.dna.save()

        self.dna.refresh_from_db()
        self.assertEqual(self.dna.reference_id, "2600000001C")
        self.assertEqual(self.dna.external_pk.code, "LIMS-EXT-1")
        self.assertIsNone(self.rna.external_pk, "A local reference works without an external one")

    def test_unnamed_extractions_are_distinct(self):
        """ unique_together is on (specimen, reference_id), and Postgres treats NULLs as distinct,
            so a specimen can hold several unnamed extractions """
        Extraction.objects.create(specimen=self.specimen)
        Extraction.objects.create(specimen=self.specimen)
        self.assertEqual(self.specimen.extraction_set.filter(reference_id__isnull=True).count(), 2)

    def test_timestamps_recorded(self):
        """ Everything in patients is a TimeStampedModel """
        for obj in (self.patient, self.specimen, self.dna):
            self.assertIsNotNone(obj.created, f"{type(obj).__name__} has no created")
            self.assertIsNotNone(obj.modified, f"{type(obj).__name__} has no modified")


class TestSpecimenTissueStatus(TestCase):
    """ variantgrid_private#2447 - a specimen nobody has described says so, rather than defaulting to
        a Germline claim that autopopulate then stamps on every classification """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("tissue_status_user", password="x")
        cls.patient = Patient.objects.create(first_name="TISSUE", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.user, cls.patient)

        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.vcf = VCF.objects.create(name="tissue_status.vcf", genotype_samples=1, genome_build=genome_build,
                                     import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())

    def test_unset_tissue_status_is_unknown(self):
        specimen = Specimen.objects.create(reference_id="UNDESCRIBED", patient=self.patient)
        self.assertEqual(specimen.tissue_status, TissueStatus.UNKNOWN)

    def _allele_origin(self, tissue_status, variants_type):
        specimen = Specimen.objects.create(reference_id=f"{tissue_status}{variants_type}",
                                           patient=self.patient, tissue_status=tissue_status)
        extraction = Extraction.objects.create(specimen=specimen, nucleic_acid_source=NucleicAcid.DNA)
        sample = Sample.objects.create(name=f"sample_{tissue_status}{variants_type}", vcf=self.vcf,
                                       patient=self.patient, extraction=extraction,
                                       variants_type=variants_type)
        data = get_evidence_fields_for_sample_and_patient(None, sample)
        return data[SpecialEKeys.ALLELE_ORIGIN]

    def test_reference_specimen_germline_call_set_is_germline(self):
        self.assertEqual(self._allele_origin(TissueStatus.REFERENCE, VariantsType.GERMLINE), "germline")

    def test_somatic_only_call_set_is_somatic_whatever_the_tissue(self):
        for tissue_status in TissueStatus:
            self.assertEqual(self._allele_origin(tissue_status, VariantsType.SOMATIC_ONLY), "somatic",
                             f"{tissue_status} + somatic only")

    def test_affected_specimen_mixed_call_set_is_unset(self):
        """ A tumour sample carrying both origins is per-variant, so the curator decides """
        self.assertIsNone(self._allele_origin(TissueStatus.AFFECTED, VariantsType.MIXED))

    def test_unknown_tissue_germline_call_set_is_unset(self):
        self.assertIsNone(self._allele_origin(TissueStatus.UNKNOWN, VariantsType.GERMLINE))

    def test_affected_specimen_germline_call_set_is_unset(self):
        self.assertIsNone(self._allele_origin(TissueStatus.AFFECTED, VariantsType.GERMLINE))

    def test_sample_without_extraction_is_unset(self):
        sample = Sample.objects.create(name="no_extraction", vcf=self.vcf, patient=self.patient,
                                       variants_type=VariantsType.GERMLINE)
        data = get_evidence_fields_for_sample_and_patient(None, sample)
        self.assertIsNone(data[SpecialEKeys.ALLELE_ORIGIN])


class TestSpecimenExtractionPermissions(TestCase):
    """ #1706 - both delegate to their patient, which is where a specimen's confidentiality lives """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.owner = User.objects.create_user("specimen_owner", password="x")
        cls.other_user = User.objects.create_user("specimen_other", password="x")
        cls.patient = Patient.objects.create(first_name="PERM", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.owner, cls.patient)

        cls.specimen = Specimen.objects.create(reference_id="PERMSPEC", patient=cls.patient)
        cls.extraction = Extraction.objects.create(specimen=cls.specimen, reference_id="PERMSPECC")

    def test_owner_sees_them(self):
        self.assertEqual(list(Specimen.filter_for_user(self.owner)), [self.specimen])
        self.assertEqual(list(Extraction.filter_for_user(self.owner)), [self.extraction])
        self.assertTrue(self.specimen.can_view(self.owner))
        self.assertTrue(self.extraction.can_view(self.owner))

    def test_other_user_sees_nothing(self):
        self.assertEqual(list(Specimen.filter_for_user(self.other_user)), [])
        self.assertEqual(list(Extraction.filter_for_user(self.other_user)), [])
        self.assertFalse(self.specimen.can_view(self.other_user))
        self.assertFalse(self.extraction.can_view(self.other_user))

    def test_external_manager_blocks_writes(self):
        """ ExternallyManagedModel.can_write is ANDed with the patient's guardian permission """
        self.assertTrue(self.specimen.can_write(self.owner))

        external_manager = ExternalModelManager.objects.create(name="read_only_lims", can_modify=False)
        self.specimen.external_pk = ExternalPK.objects.create(code="LIMS-SPEC-1", external_type="specimen",
                                                              external_manager=external_manager)
        self.specimen.save()
        self.assertFalse(self.specimen.can_write(self.owner))


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

    def test_a_second_sample_gets_its_own_extraction(self):
        """ The DNA and RNA arms of one specimen are different samples, so the second row must not
            repurpose the first sample's extraction """
        self._import_row()
        rna_sample = Sample.objects.create(name="csv_sample_rna", vcf=self.sample.vcf,
                                           import_status=ImportStatus.SUCCESS)
        assign_permission_to_user_and_groups(self.user, rna_sample)
        self._import_row(**{PatientColumns.SAMPLE_NAME: "csv_sample_rna",
                            PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE: "RNA"})

        self.assertEqual(Extraction.objects.filter(specimen__reference_id="CSVSPEC001").count(), 2)
        self.sample.refresh_from_db()
        rna_sample.refresh_from_db()
        self.assertEqual(self.sample.extraction.nucleic_acid_source, NucleicAcid.DNA)
        self.assertEqual(rna_sample.extraction.nucleic_acid_source, NucleicAcid.RNA)

    def test_blank_nucleic_acid_source_still_links_sample(self):
        self._import_row(**{PatientColumns.SPECIMEN_NUCLEIC_ACID_SOURCE: None})

        self.sample.refresh_from_db()
        self.assertEqual(self.sample.extraction.specimen.reference_id, "CSVSPEC001")
        self.assertIsNone(self.sample.extraction.nucleic_acid_source)

    def test_tissue_status_column_parses_display_value(self):
        """ parse_choice takes the key or the display value in any case """
        self._import_row(**{PatientColumns.SPECIMEN_TISSUE_STATUS: "affected / lesional"})

        specimen = Specimen.objects.get(reference_id="CSVSPEC001")
        self.assertEqual(specimen.tissue_status, TissueStatus.AFFECTED)
        record = PatientRecord.objects.get(specimen=specimen)
        self.assertEqual(record.specimen_tissue_status, TissueStatus.AFFECTED)

    def test_blank_tissue_status_column_is_unknown(self):
        self._import_row()

        specimen = Specimen.objects.get(reference_id="CSVSPEC001")
        self.assertEqual(specimen.tissue_status, TissueStatus.UNKNOWN)

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
        self.assertEqual(values[PatientColumns.SAMPLE_QUERYSET_PATH[PatientColumns.SPECIMEN_TISSUE_STATUS]],
                         TissueStatus.UNKNOWN)


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

    def test_extraction_searchable_by_its_own_reference(self):
        """ The TSO 500 container suffix lives on the extraction, not the specimen """
        self.assertIn('reference_id', ExtractionAutocompleteView.fields)

    def test_nothing_forwarded_returns_all_visible(self):
        self.assertEqual(self._queryset(SpecimenAutocompleteView).count(), 3)
        self.assertEqual(self._queryset(ExtractionAutocompleteView).count(), 4)


class TestSpecimenExtractionGrids(TestCase):
    """ The top level Specimen / Extraction grids off the patients menu (#1706) """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("grid_user", password="x")
        cls.patient = Patient.objects.create(first_name="GRID", last_name="PATIENT", patient_code="GRIDCODE")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000002", patient=cls.patient)
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000002C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        cls.rna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000002B",
                                            nucleic_acid_source=NucleicAcid.RNA)

        genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        vcf = VCF.objects.create(name="grid_test.vcf", genotype_samples=1, genome_build=genome_build,
                                 import_status=ImportStatus.SUCCESS, user=cls.user, date=timezone.now())
        cls.dna_sample = Sample.objects.create(name="grid_dna_arm", vcf=vcf, patient=cls.patient,
                                               extraction=cls.dna)
        assign_permission_to_user_and_groups(cls.user, cls.dna_sample)

    def _grid_rows(self, url_name, search=None):
        self.client.force_login(self.user)
        params = {"search[value]": search} if search else {}
        response = self.client.get(reverse(url_name), params)
        self.assertEqual(200, response.status_code)
        return {row["id"]: row for row in response.json()["data"]}

    def test_specimen_grid_counts_extractions(self):
        rows = self._grid_rows("specimen_datatables")
        self.assertEqual(rows[self.specimen.pk]["extraction_count"], 2)

    def test_extraction_grid_counts_samples_the_user_can_see(self):
        """ A sample carries its VCF's permissions rather than its patient's """
        rows = self._grid_rows("extraction_datatables")
        self.assertEqual(rows[self.dna.pk]["sample_count"], 1)
        self.assertEqual(rows[self.rna.pk]["sample_count"], 0)

        other_user = User.objects.create_user("grid_other_user", password="x")
        assign_permission_to_user_and_groups(other_user, self.patient)
        self.client.force_login(other_user)
        response = self.client.get(reverse("extraction_datatables"))
        rows = {row["id"]: row for row in response.json()["data"]}
        self.assertEqual(rows[self.dna.pk]["sample_count"], 0)

    def test_search_finds_one_arm_by_its_own_reference(self):
        rows = self._grid_rows("extraction_datatables", search="2600000002C")
        self.assertEqual(list(rows), [self.dna.pk])

    def test_search_by_patient_code(self):
        rows = self._grid_rows("specimen_datatables", search="GRIDCODE")
        self.assertEqual(list(rows), [self.specimen.pk])

    def test_unnamed_extraction_links_by_pk(self):
        """ An extraction's own reference is optional, so the link falls back to the pk """
        unnamed = Extraction.objects.create(specimen=self.specimen)
        rows = self._grid_rows("extraction_datatables")
        self.assertEqual(rows[unnamed.pk]["reference_id"]["text"], f"({unnamed.pk})")


class TestGetOrCreateExtraction(TestCase):
    """ The patient CSV names one extraction per specimen, so it has to fill in the right one without
        repurposing an arm somebody else named """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.patient = Patient.objects.create(first_name="GOC", last_name="PATIENT")
        cls.specimen = Specimen.objects.create(reference_id="2600000010", patient=cls.patient)

    def test_names_an_extraction_that_has_no_source_yet(self):
        unnamed = Extraction.objects.create(specimen=self.specimen)
        self.assertEqual(self.specimen.get_or_create_extraction(NucleicAcid.DNA), unnamed)
        unnamed.refresh_from_db()
        self.assertEqual(unnamed.nucleic_acid_source, NucleicAcid.DNA)
        self.assertEqual(self.specimen.extraction_set.count(), 1)

    def test_leaves_another_arm_alone(self):
        """ An RNA extraction added by hand must not become the CSV's DNA one """
        rna = Extraction.objects.create(specimen=self.specimen, nucleic_acid_source=NucleicAcid.RNA)
        dna = self.specimen.get_or_create_extraction(NucleicAcid.DNA)

        self.assertNotEqual(dna, rna)
        rna.refresh_from_db()
        self.assertEqual(rna.nucleic_acid_source, NucleicAcid.RNA)
        self.assertEqual(self.specimen.extraction_set.count(), 2)

    def test_no_source_reuses_whatever_is_there(self):
        """ A CSV with a blank nucleic acid column describes the extraction already on the specimen """
        dna = Extraction.objects.create(specimen=self.specimen, nucleic_acid_source=NucleicAcid.DNA)
        self.assertEqual(self.specimen.get_or_create_extraction(), dna)
        self.assertEqual(self.specimen.extraction_set.count(), 1)


class TestExtractionForm(TestCase):
    """ The specimen is fixed by the page rather than posted, so the form carries the checks Django
        skips when a unique_together field is left out """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.patient = Patient.objects.create(first_name="FORM", last_name="PATIENT")
        cls.specimen = Specimen.objects.create(reference_id="2600000020", patient=cls.patient)
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="2600000020C",
                                            nucleic_acid_source=NucleicAcid.DNA)

    def _form(self, reference_id, instance=None):
        data = {"reference_id": reference_id, "nucleic_acid_source": NucleicAcid.RNA}
        return ExtractionForm(data, instance=instance or Extraction(specimen=self.specimen))

    def test_blank_reference_saves_as_null(self):
        """ Postgres treats NULLs as distinct, which is what lets a specimen hold several unnamed
            extractions - an empty box means unnamed rather than named "" """
        form = self._form("")
        self.assertTrue(form.is_valid(), form.errors)
        extraction = form.save()
        self.assertIsNone(extraction.reference_id)

    def test_reference_already_on_this_specimen_is_rejected(self):
        form = self._form("2600000020C")
        self.assertFalse(form.is_valid())
        self.assertIn("reference_id", form.errors)

    def test_same_reference_on_another_specimen_is_fine(self):
        other = Specimen.objects.create(reference_id="2600000021", patient=self.patient)
        form = ExtractionForm({"reference_id": "2600000020C", "nucleic_acid_source": NucleicAcid.RNA},
                              instance=Extraction(specimen=other))
        self.assertTrue(form.is_valid(), form.errors)

    def test_editing_keeps_its_own_reference(self):
        form = self._form("2600000020C", instance=self.dna)
        self.assertTrue(form.is_valid(), form.errors)


class TestCreateExtractionFromSpecimenPage(TestCase):
    """ #1747 - the specimen page is where you are when you find an arm is missing """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("create_extraction_user", password="x")
        cls.patient = Patient.objects.create(first_name="CREATE", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="2600000030", patient=cls.patient)

    def _post(self, **overrides):
        self.client.force_login(self.user)
        data = {CREATE_EXTRACTION: "", "reference_id": "2600000030B",
                "nucleic_acid_source": NucleicAcid.RNA, "extraction_date": ""}
        data.update(overrides)
        return self.client.post(reverse("view_specimen", kwargs={"specimen_id": self.specimen.pk}), data)

    def test_creates_extraction_on_this_specimen(self):
        response = self._post()
        extraction = self.specimen.extraction_set.get()
        self.assertRedirects(response, extraction.get_absolute_url())
        self.assertEqual(extraction.reference_id, "2600000030B")
        self.assertEqual(extraction.nucleic_acid_source, NucleicAcid.RNA)

    def test_specimen_fields_are_left_alone(self):
        """ Both forms post to the same view, so the specimen's own form must sit out this one """
        self._post()
        self.specimen.refresh_from_db()
        self.assertEqual(self.specimen.reference_id, "2600000030")

    def test_an_externally_managed_specimen_creates_nothing(self):
        """ can_write is the external manager's answer as much as guardian's """
        manager = ExternalModelManager.objects.create(name="read_only_lims")
        self.specimen.external_pk = ExternalPK.objects.create(code="RO-1", external_type="specimen",
                                                              external_manager=manager)
        self.specimen.save()

        self._post()
        self.assertEqual(self.specimen.extraction_set.count(), 0)


class TestSpecimenExtractionSearch(TestCase):
    """ Searching a specimen reference also turns up the extractions hanging off it, which used to
        leave 2 results and nowhere to jump to """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.create_user("search_user", password="x")
        cls.patient = Patient.objects.create(first_name="SEARCH", last_name="PATIENT")
        assign_permission_to_user_and_groups(cls.user, cls.patient)
        cls.specimen = Specimen.objects.create(reference_id="SPECIMEN1", patient=cls.patient)
        # One arm built its reference off the specimen's, the other is named by the lab that took it
        cls.dna = Extraction.objects.create(specimen=cls.specimen, reference_id="SPECIMEN1C",
                                            nucleic_acid_source=NucleicAcid.DNA)
        cls.rna = Extraction.objects.create(specimen=cls.specimen, reference_id="LAB-77",
                                            nucleic_acid_source=NucleicAcid.RNA)

    def _strengths(self, search_string) -> dict:
        results = search_data(self.user, search_string, False).results
        return {r.preview.internal_url: r.match_strength for r in results}

    def test_the_specimen_reference_names_the_specimen(self):
        strengths = self._strengths("SPECIMEN1")
        self.assertEqual(strengths[self.specimen.get_absolute_url()], SearchResultMatchStrength.ID_MATCH)

    def test_an_extraction_reached_through_its_specimen_is_fuzzy(self):
        strengths = self._strengths("SPECIMEN1")
        self.assertEqual(strengths[self.rna.get_absolute_url()], SearchResultMatchStrength.FUZZY_MATCH)

    def test_an_extraction_partly_naming_itself_is_not_an_id_match(self):
        strengths = self._strengths("SPECIMEN1")
        self.assertEqual(strengths[self.dna.get_absolute_url()], SearchResultMatchStrength.STRONG_MATCH)

    def test_the_extraction_reference_names_the_extraction(self):
        strengths = self._strengths("SPECIMEN1C")
        self.assertEqual(strengths[self.dna.get_absolute_url()], SearchResultMatchStrength.ID_MATCH)

    def test_jumps_to_the_specimen(self):
        preferred = search_data(self.user, "SPECIMEN1", False).single_preferred_result()
        self.assertIsNotNone(preferred, "A whole specimen reference has somewhere to jump to")
        self.assertEqual(preferred.preview.internal_url, self.specimen.get_absolute_url())
