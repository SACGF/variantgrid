"""Import of CombinedVariantOutput.tsv - the splice VCF the loader writes from it."""
import os
import tempfile

import cyvcf2
import simplejson
from django.conf import settings
from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from genes.gene_splice import coordinate_label
from genes.models import HGNC, GeneSymbol, HGNCImport, SpliceEvent
from genes.models_enums import HGNCStatus
from library.genomics.vcf_enums import GeneIdNamespace, GeneLevelSymbolicAlt
from library.genomics.vcf_writer import percent_decode_info_value
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.models import Extraction, Patient, Specimen, SpecimenMeasure
from patients.models_enums import MatchStatus, NucleicAcid, SpecimenMeasureType
from seqauto.models import SampleFromSequencingSample, VCFFromSequencingRun
from seqauto.tests.test_extraction_link import make_sample_sheet, make_sequencing_run
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME
from snpdb.models import VCF, GenomeBuild, ImportSource, Sample, VCFSourceSettings
from snpdb.models.models_enums import ProcessingStatus
from upload.import_task_factories.import_task_factories import (
    DragenTSO500CombinedVariantOutputImportTaskFactory,
)
from upload.import_task_factories.import_task_factory import get_import_task_factories
from upload.models import (
    FileUpload,
    SimpleVCFImportInfo,
    UploadedFileTypes,
    UploadedVCF,
    UploadPipeline,
    UploadStep,
)
from upload.tasks.import_dragen_tso500_combined_variant_output_task import (
    ALT_READS_FORMAT,
    REF_READS_FORMAT,
    DragenTSO500CombinedVariantOutputCreateVCFTask,
    DragenTSO500CombinedVariantOutputInsertTask,
)
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.tso500.dragen_combined_variant_output_parser import (
    AFFECTED_EXON,
    BREAKPOINT_1,
    DNA_SAMPLE_ID,
    PAIR_ID,
    RNA_SAMPLE_ID,
    SPLICE_INFO,
    SPLICE_OBSERVATION_INFO,
    can_process_file,
    get_analysis_details,
    get_splice_rows,
    read_combined_variant_output,
)
from upload.tso500.dragen_combined_variant_output_records import (
    CombinedVariantOutputIdentityError,
    band_call,
    link_samples_to_extractions,
    link_to_sequencing_run,
    measured_date,
    parse_pair_identifiers,
    resolve_pair,
    write_specimen_measures,
)

TSO500_PAIR_DIR = os.path.join(settings.BASE_DIR, "upload", "test_data", "tso500",
                               "ExampleSample_2600000001")
COMBINED_VARIANT_OUTPUT = os.path.join(TSO500_PAIR_DIR,
                                       "ExampleSample_2600000001_CombinedVariantOutput.tsv")
EXPECTED_SPLICE_ROWS = 3


class TestCombinedVariantOutputParser(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.sections = read_combined_variant_output(COMBINED_VARIANT_OUTPUT)

    def test_analysis_details_are_key_values(self):
        details = get_analysis_details(self.sections)
        self.assertEqual("ExampleSample_RNA_2600000001B", details[RNA_SAMPLE_ID])
        self.assertEqual("2.1.1", details["Module Version"])

    def test_reads_every_splice_row(self):
        rows = get_splice_rows(self.sections)
        self.assertEqual(EXPECTED_SPLICE_ROWS, len(rows))
        self.assertEqual("chrX:66905968", rows[0][BREAKPOINT_1])

    def test_blank_affected_exon_is_none(self):
        """ AR-V7's 3' side is a cryptic exon, so the caller leaves the column empty """
        rows = {r["Gene"]: r for r in get_splice_rows(self.sections)}
        self.assertIsNone(rows["AR"][AFFECTED_EXON])
        self.assertEqual("2-7", rows["EGFR"][AFFECTED_EXON])

    def test_na_section_is_empty(self):
        """ '[Sequencing Run Details]' is a lone NA - a section with nothing in it """
        self.assertEqual([], self.sections["Sequencing Run Details"].lines)

    def test_section_it_has_no_name_for_is_still_read(self):
        """ 2.6 renames '[Exon-Level CNVs]' and adds sections of its own, so an unknown one is data """
        rows = self.sections["Exon-Level CNVs"].rows
        self.assertEqual(["BRCA1", "BRCA2"], [r["Gene"] for r in rows])

    def test_tab_padding_is_dropped(self):
        """ Every line is padded to the widest section, 11 columns """
        self.assertEqual(["Total TMB", "7.1"], self.sections["TMB"].lines[0])

    def test_claims_the_file_over_a_gene_list(self):
        factories = {type(f).__name__: f for f in get_import_task_factories()}
        combined = factories["DragenTSO500CombinedVariantOutputImportTaskFactory"]
        gene_list = factories["GeneListImportTaskFactory"]
        user = User.objects.get_or_create(username='testuser')[0]
        self.assertGreater(combined.get_processing_ability(user, COMBINED_VARIANT_OUTPUT, "tsv"),
                           gene_list.get_processing_ability(user, COMBINED_VARIANT_OUTPUT, "tsv"))

    def test_claims_a_versioned_banner(self):
        """ 2.6 writes its version into the banner line """
        with tempfile.TemporaryDirectory() as tmp_dir:
            path = os.path.join(tmp_dir, "versioned_CombinedVariantOutput.tsv")
            with open(path, "w") as f:
                f.write("DRAGEN TruSight Oncology 500 v2.6.2 Analysis Software - Combined Variant Output\t\t\n")
            self.assertTrue(can_process_file(path))

    def test_does_not_claim_other_tsvs(self):
        """ A gene table is a tsv of gene symbols - only the banner line says this is a CVO """
        self.assertFalse(can_process_file(os.path.join(settings.BASE_DIR, "seqauto", "test_data",
                                                       "reference_data", "canonical",
                                                       "fake_kit.GeneTable.tsv")))


class TestSpliceVariantVCF(TestCase):
    """ The VCF the loader writes - what the standard insert pipeline then consumes """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.grch37()
        hgnc_import = HGNCImport.objects.create()
        cls.hgnc_ids = {}
        for pk, symbol in [(644, "AR"), (3236, "EGFR"), (7029, "MET")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
            cls.hgnc_ids[symbol] = pk

    def setUp(self):
        self.user = User.objects.get_or_create(username='testuser')[0]
        self.vcf_filename = os.path.join(settings.PRIVATE_DATA_ROOT, "splice_variants.vcf")
        self.records = self._process()

    def _process(self) -> list:
        file_upload = FileUpload.objects.create(path=COMBINED_VARIANT_OUTPUT,
                                                import_source=ImportSource.COMMAND_LINE,
                                                user=self.user,
                                                name="ExampleSample_2600000001_CombinedVariantOutput.tsv",
                                                file_type=UploadedFileTypes.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT,
                                                metadata={"genome_build": "GRCh37"})
        upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        upload_step = UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                                name="Create Splice Variant VCF", sort_order=0,
                                                input_filename=COMBINED_VARIANT_OUTPUT,
                                                output_filename=self.vcf_filename)
        rows = DragenTSO500CombinedVariantOutputCreateVCFTask.process_items(upload_step)
        self.assertEqual(EXPECTED_SPLICE_ROWS, rows)
        self.reader = cyvcf2.VCF(self.vcf_filename)
        return list(self.reader)

    def test_header_names_the_rna_sample_and_the_module(self):
        """ The splice caller runs on the RNA arm, and '##source' is what VCFSourceSettings match """
        self.assertEqual(["ExampleSample_RNA_2600000001B"], self.reader.samples)
        self.assertIn("DRAGEN TSO500 CombinedVariantOutput 2.1.1", self.reader.raw_header)

    def test_records_are_gene_level_splice_alts(self):
        self.assertEqual(EXPECTED_SPLICE_ROWS, len(self.records))
        for record in self.records:
            self.assertEqual(GENE_LEVEL_CONTIG_NAME, record.CHROM)
            kind, namespace, _gene_id, label = GeneLevelSymbolicAlt.parse(record.ALT[0])
            self.assertEqual(GeneLevelSymbolicAlt.SPLICE, kind)
            self.assertEqual(GeneIdNamespace.HGNC, namespace)
            self.assertTrue(label)

    def test_seeded_junctions_get_their_label(self):
        """ The canonical label, upper-cased on the alt as a Sequence is """
        alts = {record.ALT[0] for record in self.records}
        self.assertEqual({f"<SPLICE:HGNC:{self.hgnc_ids['AR']}:V_7>",
                          f"<SPLICE:HGNC:{self.hgnc_ids['EGFR']}:V_III>",
                          f"<SPLICE:HGNC:{self.hgnc_ids['MET']}:EXON_14_SKIPPING>"}, alts)

    def _by_splice(self) -> dict:
        """ INFO values are stored as the VCF wrote them - a space is percent encoded, so decode """
        return {percent_decode_info_value(record.INFO.get(SPLICE_INFO)): record
                for record in self.records}

    def test_position_is_the_gene(self):
        self.assertEqual(self.hgnc_ids["AR"], self._by_splice()["AR-V7 splice"].POS)

    def test_read_support_is_the_junction_and_the_reference_transcript(self):
        """ A caller asserts the junction is present, so there is no GT - the sample column holds
            how many reads crossed it and how many crossed the reference transcript """
        support = {}
        for splice, record in self._by_splice().items():
            self.assertEqual([ALT_READS_FORMAT, REF_READS_FORMAT], record.FORMAT)
            support[splice] = (int(record.format(ALT_READS_FORMAT).flatten()[0]),
                               int(record.format(REF_READS_FORMAT).flatten()[0]))
        self.assertEqual({"AR-V7 splice": (27, 573), "EGFRvIII splice": (64, 1), "MET exon 14 skipping": (91, 1)}, support)

    def test_the_callers_row_rides_along_in_info(self):
        encoded = self._by_splice()["MET exon 14 skipping"].INFO.get(SPLICE_OBSERVATION_INFO)
        observation = simplejson.loads(percent_decode_info_value(encoded))
        self.assertEqual("chr7:116411708", observation[BREAKPOINT_1])
        self.assertEqual("14", observation[AFFECTED_EXON])

    def test_unnamed_junction_is_labelled_with_its_breakpoints(self):
        """ A junction no SpliceEvent names still imports - the label is its breakpoints in the
            build they were called in, which reads as raw coordinates on a report """
        splice_event = SpliceEvent.objects.get(genome_build=self.genome_build, label="v_7")
        contig = splice_event.contig
        splice_event.delete()

        records = self._process()
        expected = coordinate_label(self.genome_build, contig, 66905968, 66914514)
        self.assertIn(f"<SPLICE:HGNC:{self.hgnc_ids['AR']}:{expected.upper()}>",
                      {record.ALT[0] for record in records})

    def test_source_settings_bind_read_support_as_depth(self):
        """ The ^DRAGEN TSO500 CombinedVariantOutput row makes ALT_READS/REF_READS the depths the
            sample node and VAF use, and with no GT in the header nothing binds a genotype """
        vcf = VCF(source="DRAGEN TSO500 CombinedVariantOutput 2.1.1", genotype_samples=1)
        for vss in VCFSourceSettings.get_for_source(vcf.source):
            vss.apply_sample_field_overrides(vcf)
        self.assertEqual(ALT_READS_FORMAT, vcf.alt_depth_field)
        self.assertEqual(REF_READS_FORMAT, vcf.ref_depth_field)
        self.assertTrue(vcf.has_depth)
        self.assertFalse(vcf.has_genotype)


# The bands SA Path reports against, as its settings state them
MSI_BANDS = [(30, "MSI-High"), (10, "MSI-Low"), (0, "MSS")]
TMB_BANDS = [(10, "High"), (0, "Low")]


class TestCombinedVariantOutputRecords(TestCase):
    """ The rest of the file - the pair's patient chain, the seqauto links and the measures """

    @classmethod
    def setUpTestData(cls):
        cls.sections = read_combined_variant_output(COMBINED_VARIANT_OUTPUT)
        cls.analysis_details = get_analysis_details(cls.sections)

    def setUp(self):
        self.user = User.objects.create_user(username="cvo_records_user")
        self.identifiers = parse_pair_identifiers(self.analysis_details)

    def test_pair_names_the_patient_specimen_and_both_arms(self):
        """ The patient code is the second field of Pair ID; the accession inside each sample ID is
            the specimen and the container suffix on it that arm's extraction """
        self.assertEqual("5_C0000001_FCUP_2600000001", self.identifiers.pair_id)
        self.assertEqual("C0000001", self.identifiers.patient_code)
        self.assertEqual("2600000001", self.identifiers.specimen_reference)
        self.assertEqual("2600000001C", self.identifiers.dna.extraction_reference)
        self.assertEqual(NucleicAcid.DNA, self.identifiers.dna.nucleic_acid)
        self.assertEqual("2600000001B", self.identifiers.rna.extraction_reference)
        self.assertEqual(NucleicAcid.RNA, self.identifiers.rna.nucleic_acid)

    def test_a_bare_patient_code_pair_id_is_still_a_chain(self):
        """ The lab writes the Pair ID either way - the code alone or the whole pair sample name -
            and the accession comes off the arms' sample IDs either way """
        details = dict(self.analysis_details)
        details[PAIR_ID] = "C0000001"

        identifiers = parse_pair_identifiers(details)

        self.assertEqual("C0000001", identifiers.patient_code)
        self.assertEqual("2600000001", identifiers.specimen_reference)

    def test_pair_id_the_regex_does_not_read_is_not_a_chain(self):
        """ A code the regex cannot find names a patient we cannot identify, and the whole pair ID is
            not it - one per sequencing of the patient would leave a patient per run """
        details = dict(self.analysis_details)
        details[PAIR_ID] = "_"
        with self.assertRaises(CombinedVariantOutputIdentityError):
            parse_pair_identifiers(details)

    @override_settings(TSO500_PAIR_ID_PATIENT_CODE_REGEX=None)
    def test_with_no_regex_the_whole_pair_id_is_the_patient_code(self):
        """ A lab whose Pair ID is the patient's code, which is how Illumina documents it """
        identifiers = parse_pair_identifiers(self.analysis_details)
        self.assertEqual("5_C0000001_FCUP_2600000001", identifiers.patient_code)

    @override_settings(TSO500_PAIR_ID_PATIENT_CODE_REGEX=r"^(?P<patient_code>[^-]+)-")
    def test_another_labs_naming_is_read_by_its_own_regex(self):
        details = dict(self.analysis_details)
        details[PAIR_ID] = "C0000001-5-FCUP"
        self.assertEqual("C0000001", parse_pair_identifiers(details).patient_code)

    def test_re_sequencing_comes_back_to_the_one_patient(self):
        """ A re-sequenced patient gets a new sequencing sample ID leading the same patient code,
            so the second pair resolves to the first's patient """
        first = resolve_pair(self.identifiers, self.user)

        details = dict(self.analysis_details)
        details[PAIR_ID] = "17_C0000001_RWYN_2600000002"
        details[DNA_SAMPLE_ID] = "ExampleSample_DNA_2600000002C"
        details[RNA_SAMPLE_ID] = "ExampleSample_RNA_2600000002B"
        re_sequenced = parse_pair_identifiers(details)
        self.assertEqual("C0000001", re_sequenced.patient_code)

        second = resolve_pair(re_sequenced, self.user)
        self.assertEqual(first.patient, second.patient)
        self.assertEqual(1, Patient.objects.filter(patient_code="C0000001").count())
        self.assertNotEqual(first.specimen, second.specimen)

    def test_arms_from_different_specimens_are_not_a_pair(self):
        details = dict(self.analysis_details)
        details[DNA_SAMPLE_ID] = "ExampleSample_DNA_2600000002C"
        with self.assertRaises(CombinedVariantOutputIdentityError):
            parse_pair_identifiers(details)

    def test_chain_is_created_when_nothing_is_accessioned_yet(self):
        """ A CVO arriving first leaves a stub holding only its code for the API to fill in """
        resolved = resolve_pair(self.identifiers, self.user)
        self.assertEqual("C0000001", resolved.patient.patient_code)
        self.assertEqual("2600000001", resolved.specimen.reference_id)
        self.assertEqual(resolved.patient, resolved.specimen.patient)
        by_reference = {e.reference_id: e for e in resolved.extractions.values()}
        self.assertEqual({"2600000001C", "2600000001B"}, set(by_reference))
        self.assertEqual(NucleicAcid.DNA, by_reference["2600000001C"].nucleic_acid_source)
        self.assertEqual(NucleicAcid.RNA, by_reference["2600000001B"].nucleic_acid_source)

    def test_re_analysis_reuses_the_chain(self):
        """ A re-analysis of the pair writes the same Pair ID, so a second CVO adds nothing """
        first = resolve_pair(self.identifiers, self.user)
        second = resolve_pair(self.identifiers, self.user)
        self.assertEqual(first.patient, second.patient)
        self.assertEqual(first.specimen, second.specimen)
        self.assertEqual(1, Patient.objects.filter(patient_code="C0000001").count())
        self.assertEqual(2, Extraction.objects.filter(specimen=first.specimen).count())

    def test_specimen_held_by_another_patient_is_not_taken_over(self):
        other = Patient.objects.create(patient_code="C0009999")
        assign_permission_to_user_and_groups(self.user, other)
        Specimen.objects.create(patient=other, reference_id="2600000001")
        with self.assertRaises(CombinedVariantOutputIdentityError):
            resolve_pair(self.identifiers, self.user)

    def _make_sample(self, sample_name: str) -> Sample:
        vcf = VCF.objects.create(name=sample_name, date=timezone.now(), user=self.user,
                                 genotype_samples=1, genome_build=GenomeBuild.grch37())
        assign_permission_to_user_and_groups(self.user, vcf)
        return Sample.objects.create(vcf=vcf, name=sample_name, vcf_sample_name=sample_name)

    def test_both_arms_samples_link_to_their_extraction(self):
        """ The file's sample IDs are the arms' VCF sample names, which is the whole join """
        samples = {arm.sample_id: self._make_sample(arm.sample_id) for arm in self.identifiers.arms}
        resolved = resolve_pair(self.identifiers, self.user)
        self.assertEqual(2, link_samples_to_extractions(resolved, self.user))

        for sample_id, sample in samples.items():
            sample.refresh_from_db()
            self.assertEqual(resolved.extractions[sample_id], sample.extraction)
            self.assertEqual(MatchStatus.MATCHED, sample.extraction_match_status)

    def test_the_rna_arm_links_the_vcf_to_its_sequencing_run(self):
        rna_sample_id = self.identifiers.rna.sample_id
        sequencing_run = make_sequencing_run("TSO500_CVO")
        _sample_sheet, sequencing_samples = make_sample_sheet(sequencing_run, [rna_sample_id])
        sample = self._make_sample(rna_sample_id)

        self.assertEqual(sequencing_run, link_to_sequencing_run(sample.vcf, sample, rna_sample_id))
        self.assertTrue(VCFFromSequencingRun.objects.filter(vcf=sample.vcf,
                                                            sequencing_run=sequencing_run).exists())
        self.assertTrue(SampleFromSequencingSample.objects.filter(
            sample=sample, sequencing_sample=sequencing_samples[0]).exists())

    def test_a_sample_sheet_without_the_arm_leaves_no_link_rows(self):
        sample = self._make_sample(self.identifiers.rna.sample_id)
        make_sample_sheet(make_sequencing_run("TSO500_OTHER"), ["SOMETHING_ELSE"])
        self.assertIsNone(link_to_sequencing_run(sample.vcf, sample, self.identifiers.rna.sample_id))
        self.assertFalse(VCFFromSequencingRun.objects.filter(vcf=sample.vcf).exists())

    def test_the_pairs_measures_are_written_against_specimen_and_dna_arm(self):
        resolved = resolve_pair(self.identifiers, self.user)
        measures = write_specimen_measures(self.sections, resolved, self.identifiers, self.user,
                                           method="DRAGEN TSO500 CombinedVariantOutput 2.1.1",
                                           date=measured_date(self.analysis_details))
        by_type = {m.measure_type: m for m in measures}
        self.assertEqual({SpecimenMeasureType.TMB: 7.1, SpecimenMeasureType.MSI: 2.48,
                          SpecimenMeasureType.GIS: 31.0, SpecimenMeasureType.TUMOUR_FRACTION: 0.62,
                          SpecimenMeasureType.PLOIDY: 2.10},
                         {measure_type: m.value for measure_type, m in by_type.items()})
        tmb = by_type[SpecimenMeasureType.TMB]
        self.assertEqual(resolved.specimen, tmb.specimen)
        self.assertEqual(resolved.extractions[self.identifiers.dna.sample_id], tmb.extraction)
        self.assertEqual("mut/Mb", tmb.unit)
        self.assertEqual("1.27", tmb.source_payload["Coding Region Size in Megabases"])

    def _write_measures(self) -> dict:
        resolved = resolve_pair(self.identifiers, self.user)
        measures = write_specimen_measures(self.sections, resolved, self.identifiers, self.user,
                                           method="DRAGEN TSO500 CombinedVariantOutput 2.1.1")
        return {m.measure_type: m for m in measures}

    @override_settings(TSO500_MSI_MIN_USABLE_SITES=40, TSO500_MSI_CALL_BANDS=MSI_BANDS,
                       TSO500_TMB_CALL_BANDS=TMB_BANDS)
    def test_msi_and_tmb_are_called_against_the_labs_bands(self):
        """ 121 usable sites is enough to call MSI, and 2.48% unstable is MSS; 7.1 mut/Mb is Low """
        by_type = self._write_measures()

        msi = by_type[SpecimenMeasureType.MSI]
        self.assertEqual("MSS", msi.call)
        self.assertEqual("MSI-High >= 30%, MSI-Low >= 10%, MSS < 10% unstable sites, needs >= 40 usable sites",
                         msi.threshold)
        self.assertEqual("Low", by_type[SpecimenMeasureType.TMB].call)
        self.assertEqual("High >= 10 mut/Mb, Low < 10 mut/Mb", by_type[SpecimenMeasureType.TMB].threshold)
        # The measures the lab has no policy for are the number alone
        self.assertIsNone(by_type[SpecimenMeasureType.GIS].call)

    def test_band_call_is_the_first_lower_bound_reached(self):
        self.assertEqual("MSI-High", band_call(30, MSI_BANDS))
        self.assertEqual("MSI-Low", band_call(29.9, MSI_BANDS))
        self.assertEqual("MSS", band_call(0, MSI_BANDS))
        self.assertEqual("High", band_call(10, TMB_BANDS))

    @override_settings(TSO500_MSI_MIN_USABLE_SITES=200, TSO500_MSI_CALL_BANDS=MSI_BANDS)
    def test_too_few_usable_msi_sites_cannot_be_called(self):
        """ The percentage means nothing off 121 sites when the lab wants 200 - the threshold that
            was applied is still recorded, so 'why is there no call' stays answerable """
        msi = self._write_measures()[SpecimenMeasureType.MSI]

        self.assertIsNone(msi.call)
        self.assertIn("needs >= 200 usable sites", msi.threshold)

    def test_thresholds_unset_writes_the_value_and_no_call(self):
        """ The policy is the lab's - an installation without one sends the measure as it always did """
        by_type = self._write_measures()

        self.assertIsNone(by_type[SpecimenMeasureType.MSI].call)
        self.assertIsNone(by_type[SpecimenMeasureType.TMB].call)
        self.assertIsNone(by_type[SpecimenMeasureType.TMB].threshold)

    def test_re_analysis_replaces_the_measures(self):
        """ One current value per measure - the report wants a single TMB, not a history """
        resolved = resolve_pair(self.identifiers, self.user)
        for _ in range(2):
            write_specimen_measures(self.sections, resolved, self.identifiers, self.user,
                                    method="DRAGEN TSO500 CombinedVariantOutput 2.1.1")
        self.assertEqual(5, SpecimenMeasure.objects.filter(specimen=resolved.specimen).count())


class TestCombinedVariantOutputInsertTask(TestCase):
    """ The step that runs once the Sample exists - the whole of the rest of the file """

    def setUp(self):
        self.user = User.objects.create_user(username="cvo_insert_user")
        self.rna_sample_id = "ExampleSample_RNA_2600000001B"
        vcf = VCF.objects.create(name=self.rna_sample_id, date=timezone.now(), user=self.user,
                                 genotype_samples=1, genome_build=GenomeBuild.grch37())
        assign_permission_to_user_and_groups(self.user, vcf)
        self.sample = Sample.objects.create(vcf=vcf, name=self.rna_sample_id,
                                            vcf_sample_name=self.rna_sample_id)
        file_upload = FileUpload.objects.create(path=COMBINED_VARIANT_OUTPUT,
                                                import_source=ImportSource.COMMAND_LINE,
                                                user=self.user,
                                                name="ExampleSample_2600000001_CombinedVariantOutput.tsv",
                                                file_type=UploadedFileTypes.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT)
        self.upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
        UploadedVCF.objects.create(file_upload=file_upload, upload_pipeline=self.upload_pipeline,
                                   vcf=vcf)
        self.upload_step = UploadStep.objects.create(upload_pipeline=self.upload_pipeline,
                                                     name="CombinedVariantOutput records",
                                                     sort_order=1)

    def test_the_pair_lands_with_its_sample_and_measures(self):
        measures = DragenTSO500CombinedVariantOutputInsertTask.process_items(self.upload_step)
        self.assertEqual(5, measures)

        patient = Patient.objects.get(patient_code="C0000001")
        specimen = Specimen.objects.get(patient=patient, reference_id="2600000001")
        self.sample.refresh_from_db()
        self.assertEqual("2600000001B", self.sample.extraction.reference_id)
        self.assertEqual(specimen, self.sample.extraction.specimen)
        self.assertEqual(5, SpecimenMeasure.objects.filter(specimen=specimen).count())

    def test_an_unregistered_sequencing_sample_is_a_message_not_a_failure(self):
        DragenTSO500CombinedVariantOutputInsertTask.process_items(self.upload_step)
        messages = SimpleVCFImportInfo.objects.filter(upload_step__upload_pipeline=self.upload_pipeline)
        self.assertIn(self.rna_sample_id, "\n".join(m.message for m in messages))
        self.assertFalse(VCFFromSequencingRun.objects.exists())

    def test_the_seqauto_rows_are_written_for_the_rna_arm(self):
        sequencing_run = make_sequencing_run("TSO500_INSERT")
        _sample_sheet, sequencing_samples = make_sample_sheet(sequencing_run, [self.rna_sample_id])

        DragenTSO500CombinedVariantOutputInsertTask.process_items(self.upload_step)
        self.assertTrue(VCFFromSequencingRun.objects.filter(vcf=self.sample.vcf,
                                                            sequencing_run=sequencing_run).exists())
        self.assertTrue(SampleFromSequencingSample.objects.filter(
            sample=self.sample, sequencing_sample=sequencing_samples[0]).exists())

    def test_a_pair_with_no_splice_calls_still_gets_its_records(self):
        """ Most pairs have no splice call, which is an empty VCF - and that skips every step
            waiting on data insertion """
        self.upload_step.delete()
        DragenTSO500CombinedVariantOutputImportTaskFactory().create_import_task(self.upload_pipeline)
        steps = {step.name: step for step in self.upload_pipeline.uploadstep_set.all()}
        ImportVCFStepTask._handle_no_vcf_records(steps["Schedule Parallel VCF Processing tasks"])

        statuses = dict(self.upload_pipeline.uploadstep_set.values_list("name", "status"))
        self.assertEqual(ProcessingStatus.SKIPPED, statuses["VCFCheckAnnotationTask"])
        self.assertEqual(ProcessingStatus.CREATED, statuses[DragenTSO500CombinedVariantOutputInsertTask.__name__])
