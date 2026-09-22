"""Import of MetricsOutput.tsv - the LibraryQC rows the loader writes from it (sapath#455)."""
import os
import tempfile

from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings

from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Patient, Specimen
from patients.models_enums import MatchStatus, NucleicAcid
from patients.tasks.extraction_matching_tasks import reconcile_pending_extractions
from seqauto.models import (
    SAMPLE_SHEET_PAIR_ID_COLUMN,
    SAMPLE_SHEET_SAMPLE_TYPE_COLUMN,
    LibraryQC,
    SequencingSampleData,
)
from seqauto.models.models_enums import LibraryQCCategory
from seqauto.tests.test_extraction_link import make_sample_sheet, make_sequencing_run
from snpdb.models import ImportSource
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import FileUpload, UploadedFileTypes, UploadPipeline
from upload.tasks.import_dragen_tso500_metrics_output_task import (
    ImportDragenTSO500MetricsOutputTask,
)
from upload.tso500.dragen_metrics_output_parser import (
    GUIDELINE_SETTING,
    can_process_file,
    read_library_qc,
    read_sections,
)
from upload.tso500.dragen_metrics_output_records import (
    metrics_measured_date,
    metrics_method,
    write_library_qc,
)

TSO500_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500")
TSO500_PAIR_DIR = os.path.join(TSO500_DIR, "ExampleSample_2600000001")
# DRAGEN 2.1.1's per-pair copy, whose columns are the two arms' sample IDs
PAIR_METRICS_OUTPUT = os.path.join(TSO500_PAIR_DIR, "ExampleSample_2600000001_MetricsOutput.tsv")
# The run-level 2.6.2 file the pipeline sends - a column per pair, carrying both of its arms
RUN_METRICS_OUTPUT = os.path.join(TSO500_DIR, "MetricsOutput_orig.tsv")
SEQUENCING_RUN_NAME = "260101_M02027_0001_000000000-TSO500"
DNA_SAMPLE_ID = "ExampleSample_DNA_2600000001C"
RNA_SAMPLE_ID = "ExampleSample_RNA_2600000001B"
FIRST_PAIR = "5_C0000001_FCUP_2600000001"
SECOND_PAIR = "7_C0000002_ABCD_2600000002"
CONTROL_PAIR = "9_0PRI_2600000003"
# The lab writes the Pair ID either way - this one is the patient code on its own, with no accession
BARE_PAIR = "C0000004"
BARE_PAIR_SHEET_SAMPLE = "1_TSO_DNAHRD_C0000004_2600000004C_B4"


def make_pair_sample_sheet(sequencing_run, arms: dict[str, dict[str, str]], sheet_hash="HASH1"):
    """ A TSO500 sheet as the pipeline posts it: each arm's row carries the sheet's Pair_ID and
        Sample_Type columns as SequencingSampleData. arms is {pair_id: {'DNA': name, 'RNA': name}} """
    names = [name for pair_arms in arms.values() for name in pair_arms.values()]
    _, sequencing_samples = make_sample_sheet(sequencing_run, names, sheet_hash=sheet_hash)
    by_name = {ss.sample_name: ss for ss in sequencing_samples}
    for pair_id, pair_arms in arms.items():
        for sample_type, name in pair_arms.items():
            SequencingSampleData.objects.create(sequencing_sample=by_name[name],
                                                column=SAMPLE_SHEET_PAIR_ID_COLUMN, value=pair_id)
            SequencingSampleData.objects.create(sequencing_sample=by_name[name],
                                                column=SAMPLE_SHEET_SAMPLE_TYPE_COLUMN, value=sample_type)
    return by_name
NO_ACCESSION_PAIR = "11_NTC"
ALL_PAIRS = (FIRST_PAIR, SECOND_PAIR, CONTROL_PAIR, BARE_PAIR, NO_ACCESSION_PAIR)
DNA_CATEGORIES = {LibraryQCCategory.DNA, LibraryQCCategory.SMALL_VARIANT_TMB,
                  LibraryQCCategory.MSI, LibraryQCCategory.CNV, LibraryQCCategory.GIS}
ALL_CATEGORIES = DNA_CATEGORIES | {LibraryQCCategory.RNA}


def edited_metrics_file(tmp_dir: str, *replacements) -> str:
    """ The run-level fixture with a line or two changed - a lab file whose numbers we control """
    with open(RUN_METRICS_OUTPUT, encoding="utf-8-sig") as f:
        contents = f.read()
    for old, new in replacements:
        assert old in contents, old
        contents = contents.replace(old, new)
    path = os.path.join(tmp_dir, "edited_MetricsOutput.tsv")
    with open(path, "w") as f:
        f.write(contents)
    return path


def read_pairs(path: str) -> dict:
    return {library.pair_id: library for library in read_library_qc(read_sections(path))}


class TestMetricsOutputParser(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.pairs = read_pairs(RUN_METRICS_OUTPUT)

    def test_a_column_carries_both_arms_of_its_pair(self):
        """ 2.6.2's column is the CVO's Pair ID, so the DNA and the RNA sections both have values in it """
        first = self.pairs[FIRST_PAIR]

        self.assertEqual(dict.fromkeys(ALL_CATEGORIES, True),
                         {category: first.category_passed(category) for category in ALL_CATEGORIES})

    def test_a_metric_outside_its_guideline_fails_its_category(self):
        second = self.pairs[SECOND_PAIR]

        self.assertFalse(second.category_passed(LibraryQCCategory.SMALL_VARIANT_TMB))
        self.assertFalse(second.categories[LibraryQCCategory.SMALL_VARIANT_TMB]["MEDIAN_EXON_COVERAGE"]["passed"])
        self.assertTrue(second.category_passed(LibraryQCCategory.MSI))

    def test_an_na_guideline_is_no_bound(self):
        """ USABLE_MSI_SITES is LSL 40 and USL NA, so there is no such thing as too many """
        metric = self.pairs[FIRST_PAIR].categories[LibraryQCCategory.MSI]["USABLE_MSI_SITES"]

        self.assertEqual(40, metric["lsl"])
        self.assertIsNone(metric["usl"])
        self.assertTrue(metric["passed"])

    def test_a_metric_the_version_added_counts_towards_its_category(self):
        """ 2.6.2 added these two, and only the section names are keyed on here """
        first = self.pairs[FIRST_PAIR]

        self.assertIn("PCT_CHIMERIC_READS", first.categories[LibraryQCCategory.SMALL_VARIANT_TMB])
        self.assertIn("EXCESSIVE_TF", first.categories[LibraryQCCategory.GIS])

    def test_the_metric_column_carries_its_unit(self):
        metrics = self.pairs[FIRST_PAIR].categories[LibraryQCCategory.SMALL_VARIANT_TMB]

        self.assertEqual("bp", metrics["MEDIAN_INSERT_SIZE"]["unit"])
        self.assertEqual("%", metrics["PCT_EXON_50X"]["unit"])
        # '(NA)' is a metric with no unit at all
        self.assertIsNone(self.pairs[FIRST_PAIR].categories[LibraryQCCategory.DNA]["CONTAMINATION_SCORE"]["unit"])

    def test_a_run_that_did_not_finish_for_a_pair(self):
        self.assertTrue(self.pairs[FIRST_PAIR].completed)
        self.assertFalse(self.pairs[SECOND_PAIR].completed)

    @override_settings(TSO500_LIBRARY_QC_GUIDELINES={
        ("RNA Library QC Metrics", "TOTAL_ON_TARGET_READS"): (20000000, None)})
    def test_the_labs_own_guideline_overrides_the_files(self):
        """ The 2.6.2 file's LSL is 2.5M where a lab's methods paragraph may quote its own number """
        metric = read_pairs(RUN_METRICS_OUTPUT)[FIRST_PAIR].categories[LibraryQCCategory.RNA]["TOTAL_ON_TARGET_READS"]

        self.assertEqual(20000000, metric["lsl"])
        self.assertFalse(metric["passed"])
        self.assertEqual(GUIDELINE_SETTING, metric["guideline_source"])
        # MEDIAN_INSERT_SIZE is in two sections, so an override is keyed on the section too
        self.assertEqual("file", read_pairs(RUN_METRICS_OUTPUT)[FIRST_PAIR]
                         .categories[LibraryQCCategory.SMALL_VARIANT_TMB]["MEDIAN_INSERT_SIZE"]["guideline_source"])

    def test_the_older_per_pair_files_sample_columns_read_the_same_way(self):
        """ 2.1.1 wrote a column per sample rather than per pair - a column is a column """
        libraries = read_pairs(PAIR_METRICS_OUTPUT)

        self.assertTrue(libraries[DNA_SAMPLE_ID].category_passed(LibraryQCCategory.CNV))
        self.assertIsNone(libraries[DNA_SAMPLE_ID].category_passed(LibraryQCCategory.RNA))
        self.assertTrue(libraries[RNA_SAMPLE_ID].category_passed(LibraryQCCategory.RNA))

    def test_claims_the_file_over_a_gene_list_and_the_combined_variant_output(self):
        factories = {type(f).__name__: f for f in get_import_task_factories()}
        user = User.objects.get_or_create(username='testuser')[0]
        abilities = {name: factories[name].get_processing_ability(user, RUN_METRICS_OUTPUT, "tsv")
                     for name in ("DragenTSO500MetricsOutputImportTaskFactory",
                                  "DragenTSO500CombinedVariantOutputImportTaskFactory",
                                  "GeneListImportTaskFactory")}

        self.assertGreater(abilities["DragenTSO500MetricsOutputImportTaskFactory"],
                           abilities["GeneListImportTaskFactory"])
        self.assertEqual(0, abilities["DragenTSO500CombinedVariantOutputImportTaskFactory"])

    def test_does_not_claim_the_combined_variant_output(self):
        """ The two files share a layout - only the banner line says which is which """
        combined = os.path.join(TSO500_PAIR_DIR, "ExampleSample_2600000001_CombinedVariantOutput.tsv")
        self.assertFalse(can_process_file(combined))


class TestMetricsOutputRecords(TestCase):
    """ The rows a run's file writes, and the specimen each pair claims by its accession """

    @classmethod
    def setUpTestData(cls):
        cls.sections = read_sections(RUN_METRICS_OUTPUT)
        cls.libraries = read_library_qc(cls.sections)

    def setUp(self):
        self.user = User.objects.create_user(username="metrics_records_user")
        self.patient = Patient.objects.create(patient_code="C0000001")
        assign_permission_to_user_and_groups(self.user, self.patient)

    def _specimen(self, reference_id="2600000001") -> Specimen:
        return Specimen.objects.create(patient=self.patient, reference_id=reference_id)

    def _write(self, sequencing_run_name=SEQUENCING_RUN_NAME) -> dict:
        rows = write_library_qc(self.libraries, self.user, sequencing_run_name,
                                method=metrics_method(self.sections),
                                date=metrics_measured_date(self.sections))
        return {(row.pair_id, row.category): row for row in rows}

    def test_one_row_per_pair_per_category_it_has_qc_for(self):
        rows = self._write()

        self.assertEqual({(pair, category) for pair in ALL_PAIRS for category in ALL_CATEGORIES},
                         set(rows))
        cnv = rows[(FIRST_PAIR, LibraryQCCategory.CNV)]
        self.assertTrue(cnv.passed)
        self.assertEqual(NucleicAcid.DNA, cnv.nucleic_acid)
        self.assertEqual(NucleicAcid.RNA, rows[(FIRST_PAIR, LibraryQCCategory.RNA)].nucleic_acid)
        self.assertEqual("DRAGEN TSO500 MetricsOutput 2.6.2.4", cnv.method)
        # The whole section, so which number the call came off stays answerable
        self.assertIn("GENE_SCALED_MAD", cnv.metrics)

    def test_a_column_nothing_on_the_run_is_named_for_is_parked(self):
        """ A run holds more than its TSO 500 pairs - a control's QC is kept, attached to nothing """
        row = self._write()[(NO_ACCESSION_PAIR, LibraryQCCategory.CNV)]

        self.assertEqual("", row.specimen_reference)
        self.assertEqual(MatchStatus.PENDING, row.specimen_match_status)
        self.assertIn(SEQUENCING_RUN_NAME, row.specimen_match_error)

    def test_a_bare_patient_code_column_takes_its_accession_off_the_run(self):
        """ The sheet's sample names are the one place the code and the accession are written together """
        sequencing_run = make_sequencing_run(SEQUENCING_RUN_NAME)
        make_sample_sheet(sequencing_run, [BARE_PAIR_SHEET_SAMPLE])
        specimen = self._specimen("2600000004")

        row = self._write()[(BARE_PAIR, LibraryQCCategory.CNV)]

        self.assertEqual("2600000004", row.specimen_reference)
        self.assertEqual(specimen, row.specimen)
        self.assertEqual(MatchStatus.MATCHED, row.specimen_match_status)

    def test_a_bare_patient_code_column_with_no_sample_sheet_is_parked(self):
        row = self._write()[(BARE_PAIR, LibraryQCCategory.CNV)]

        self.assertEqual("", row.specimen_reference)
        self.assertEqual(MatchStatus.PENDING, row.specimen_match_status)
        self.assertIn(BARE_PAIR, row.specimen_match_error)

    def test_the_second_pairs_failing_metric_and_unfinished_run(self):
        rows = self._write()

        variants = rows[(SECOND_PAIR, LibraryQCCategory.SMALL_VARIANT_TMB)]
        self.assertFalse(variants.passed)
        self.assertFalse(variants.metrics["MEDIAN_EXON_COVERAGE"]["passed"])
        self.assertTrue(rows[(SECOND_PAIR, LibraryQCCategory.MSI)].passed)
        self.assertFalse(rows[(SECOND_PAIR, LibraryQCCategory.RNA)].completed)
        self.assertEqual("run did not complete",
                         rows[(SECOND_PAIR, LibraryQCCategory.RNA)].status_description)

    def test_a_pair_whose_specimen_exists_is_matched(self):
        specimen = self._specimen()

        row = self._write()[(FIRST_PAIR, LibraryQCCategory.CNV)]

        self.assertEqual(specimen, row.specimen)
        self.assertEqual(MatchStatus.MATCHED, row.specimen_match_status)

    def test_a_pair_with_no_specimen_yet_is_parked_and_reconciled(self):
        """ Accessioning a case is the CombinedVariantOutput's job, so QC arriving first parks """
        row = self._write()[(FIRST_PAIR, LibraryQCCategory.CNV)]
        self.assertIsNone(row.specimen)
        self.assertEqual(MatchStatus.PENDING, row.specimen_match_status)
        self.assertEqual("2600000001", row.specimen_reference)

        specimen = self._specimen()
        reconcile_pending_extractions()

        row.refresh_from_db()
        self.assertEqual(specimen, row.specimen)
        self.assertEqual(MatchStatus.MATCHED, row.specimen_match_status)

    def test_a_pair_with_no_patient_code_claims_its_specimen_like_any_other(self):
        row = self._write()[(CONTROL_PAIR, LibraryQCCategory.CNV)]

        self.assertEqual("2600000003", row.specimen_reference)
        self.assertEqual(MatchStatus.PENDING, row.specimen_match_status)

    def test_the_same_pair_on_another_run_is_a_second_row(self):
        """ A pair ID's leading sequencing number changes on re-sequencing, but nothing promises it does """
        self._write()
        self._write("260201_M02027_0002_000000000-TSO500")

        self.assertEqual(2, LibraryQC.objects.filter(pair_id=FIRST_PAIR,
                                                     category=LibraryQCCategory.CNV).count())

    def test_a_re_analysis_of_the_run_replaces_its_rows(self):
        self._write()
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = edited_metrics_file(tmp_dir, ("GENE_SCALED_MAD (Count)\t0\t0.134\t0.047",
                                                 "GENE_SCALED_MAD (Count)\t0\t0.134\t0.9"))
            sections = read_sections(path)
            write_library_qc(read_library_qc(sections), self.user, SEQUENCING_RUN_NAME,
                             method=metrics_method(sections))

        self.assertEqual(30, LibraryQC.objects.count())
        cnv = LibraryQC.objects.get(pair_id=FIRST_PAIR, category=LibraryQCCategory.CNV)
        self.assertFalse(cnv.passed)

    def test_a_registered_run_is_linked(self):
        sequencing_run = make_sequencing_run(SEQUENCING_RUN_NAME)

        row = self._write()[(FIRST_PAIR, LibraryQCCategory.CNV)]

        self.assertEqual(sequencing_run, row.sequencing_run)

    def test_a_run_seqauto_has_not_seen_leaves_the_link_null(self):
        """ The rows are still written and still keyed - the run can be registered afterwards """
        row = self._write()[(FIRST_PAIR, LibraryQCCategory.CNV)]

        self.assertIsNone(row.sequencing_run)
        self.assertIsNone(row.sequencing_sample)
        self.assertEqual(SEQUENCING_RUN_NAME, row.sequencing_run_name)

    def test_each_arm_links_its_row_on_the_sheet(self):
        """ The sheet's Pair_ID / Sample_Type columns say which sample is which arm of the pair """
        sequencing_run = make_sequencing_run(SEQUENCING_RUN_NAME)
        by_name = make_pair_sample_sheet(sequencing_run, {FIRST_PAIR: {"DNA": DNA_SAMPLE_ID, "RNA": RNA_SAMPLE_ID}})

        rows = self._write()

        self.assertEqual(by_name[DNA_SAMPLE_ID], rows[(FIRST_PAIR, LibraryQCCategory.CNV)].sequencing_sample)
        self.assertEqual(by_name[DNA_SAMPLE_ID], rows[(FIRST_PAIR, LibraryQCCategory.MSI)].sequencing_sample)
        self.assertEqual(by_name[RNA_SAMPLE_ID], rows[(FIRST_PAIR, LibraryQCCategory.RNA)].sequencing_sample)
        self.assertIsNone(rows[(SECOND_PAIR, LibraryQCCategory.CNV)].sequencing_sample)

    def test_a_bare_patient_code_column_takes_its_accession_off_its_linked_arm(self):
        sequencing_run = make_sequencing_run(SEQUENCING_RUN_NAME)
        by_name = make_pair_sample_sheet(sequencing_run, {BARE_PAIR: {"DNA": BARE_PAIR_SHEET_SAMPLE}})
        specimen = self._specimen("2600000004")

        row = self._write()[(BARE_PAIR, LibraryQCCategory.CNV)]

        self.assertEqual(by_name[BARE_PAIR_SHEET_SAMPLE], row.sequencing_sample)
        self.assertEqual(specimen, row.specimen)

    def test_a_sheet_registered_after_the_file_is_linked_by_reconcile(self):
        row = self._write()[(FIRST_PAIR, LibraryQCCategory.RNA)]
        self.assertIsNone(row.sequencing_run)

        sequencing_run = make_sequencing_run(SEQUENCING_RUN_NAME)
        by_name = make_pair_sample_sheet(sequencing_run, {FIRST_PAIR: {"DNA": DNA_SAMPLE_ID, "RNA": RNA_SAMPLE_ID}})
        reconcile_pending_extractions()

        row.refresh_from_db()
        self.assertEqual(sequencing_run, row.sequencing_run)
        self.assertEqual(by_name[RNA_SAMPLE_ID], row.sequencing_sample)

    def test_a_re_sent_sheet_moves_the_link_to_its_new_rows(self):
        sequencing_run = make_sequencing_run(SEQUENCING_RUN_NAME)
        make_pair_sample_sheet(sequencing_run, {FIRST_PAIR: {"DNA": DNA_SAMPLE_ID}})
        row = self._write()[(FIRST_PAIR, LibraryQCCategory.CNV)]

        by_name = make_pair_sample_sheet(sequencing_run, {FIRST_PAIR: {"DNA": DNA_SAMPLE_ID}}, sheet_hash="HASH2")
        reconcile_pending_extractions()

        row.refresh_from_db()
        self.assertEqual(by_name[DNA_SAMPLE_ID], row.sequencing_sample)


class TestMetricsOutputImportTask(TestCase):
    """ The single shot import - the file has no variants, so there is no VCF pipeline """

    def setUp(self):
        self.user = User.objects.create_user(username="metrics_import_user")
        patient = Patient.objects.create(patient_code="C0000001")
        assign_permission_to_user_and_groups(self.user, patient)
        self.specimen = Specimen.objects.create(patient=patient, reference_id="2600000001")

    def _file_upload(self, metadata) -> FileUpload:
        file_upload = FileUpload.objects.create(path=RUN_METRICS_OUTPUT,
                                                import_source=ImportSource.COMMAND_LINE,
                                                user=self.user, name="MetricsOutput_orig.tsv",
                                                file_type=UploadedFileTypes.DRAGEN_TSO500_METRICS_OUTPUT,
                                                metadata=metadata)
        UploadPipeline.objects.create(file_upload=file_upload)
        return file_upload

    def test_the_runs_rows_land_with_the_upload_that_wrote_them(self):
        file_upload = self._file_upload({"sequencing_run": SEQUENCING_RUN_NAME})

        self.assertEqual(30, ImportDragenTSO500MetricsOutputTask.process_items(file_upload))

        self.assertEqual(30, LibraryQC.objects.filter(file_upload=file_upload).count())
        self.assertEqual(6, LibraryQC.objects.filter(specimen=self.specimen).count())
        self.assertEqual(self.user, LibraryQC.objects.first().user)

    def test_an_upload_that_names_no_run_fails_saying_so(self):
        file_upload = self._file_upload({})

        with self.assertRaises(ValueError) as cm:
            ImportDragenTSO500MetricsOutputTask.process_items(file_upload)

        self.assertIn("sequencing_run", str(cm.exception))
        self.assertFalse(LibraryQC.objects.exists())
