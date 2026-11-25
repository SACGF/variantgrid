"""
#2898 - auto-launch rules matched on variant caller / VCF source, and launched per extraction / specimen /
patient once every file registered for the group on a run has imported
"""
from unittest.mock import patch

from django.contrib.auth.models import User
from django.test import TestCase, override_settings
from django.utils import timezone

from analysis.auto_launch import (
    AutoLaunchSubject,
    auto_launch_group_analyses_for_sample,
    get_auto_launch_analysis_template_matches,
)
from analysis.enums import AnalysisTemplateRunTrigger, AnalysisTemplateType
from analysis.forms.forms import AutoLaunchAnalysisTemplateForm
from analysis.models.models_analysis import (
    Analysis,
    AnalysisTemplate,
    AnalysisVariable,
    AutoLaunchAnalysisTemplate,
    GroupAnalysisTemplateRun,
)
from analysis.models.nodes.sources.sample_node import SampleNode
from annotation.fake_data import get_fake_annotation_version
from eventlog.models import Event
from library.guardian_utils import assign_permission_to_user_and_groups
from patients.enums import NucleicAcid, SampleSourceLevel
from patients.models import Extraction, Patient, Specimen
from seqauto.enums import DataGeneration
from seqauto.models import (
    SampleSheet,
    Sequencer,
    SequencerModel,
    SequencingRun,
    SequencingRunCurrentSampleSheet,
    SequencingSample,
    SingleSampleVCF,
    VariantCaller,
    VCFFromSequencingRun,
)
from seqauto.models.models_upload import BackendVCF
from snpdb.enums import ImportStatus
from snpdb.models import VCF, Cohort, CohortGenotypeCollection, CohortSample, GenomeBuild, Sample
from upload.models import FileUpload, UploadedVCF

DNA_CALLERS = ("TSO_dragen_caller", "TSO_dragen_cnv_caller")
RNA_CALLERS = ("TSO_fusion_caller",)


@override_settings(ANALYSIS_NODE_CACHE_Q=False)
class AutoLaunchTestBase(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="test_auto_launch_user")[0]
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.callers = {name: VariantCaller.objects.create(name=name, version="2.6")
                       for name in DNA_CALLERS + RNA_CALLERS}
        cls.template = cls._create_template("auto launch specimen template", SampleSourceLevel.SPECIMEN)
        AutoLaunchAnalysisTemplate.objects.create(template=cls.template, launch_level=SampleSourceLevel.SPECIMEN)

    @classmethod
    def _create_template(cls, name, level) -> AnalysisTemplate:
        analysis = Analysis(genome_build=cls.grch37, template_type=AnalysisTemplateType.TEMPLATE)
        analysis.set_defaults_and_save(cls.user)
        field = SampleNode.SOURCE_LEVEL_FIELDS[level]
        node = SampleNode.objects.create(analysis=analysis, source_level=level)
        AnalysisVariable.objects.create(node=node, field=field, class_name=f"patients.{field.capitalize()}")
        template = AnalysisTemplate.objects.create(name=name, user=cls.user, analysis=analysis)
        AnalysisTemplate.objects.get(pk=template.pk).new_version("%(template)s for %(input)s").activate()
        return template

    def setUp(self):
        patient = Patient.objects.create(first_name="AUTO", last_name="LAUNCH")
        assign_permission_to_user_and_groups(self.user, patient)
        self.specimen = Specimen.objects.create(patient=patient, reference_id="2600000001")
        self.dna = Extraction.objects.create(specimen=self.specimen, reference_id="2600000001C",
                                             nucleic_acid_source=NucleicAcid.DNA)
        self.rna = Extraction.objects.create(specimen=self.specimen, reference_id="2600000001D",
                                             nucleic_acid_source=NucleicAcid.RNA)

    def _make_run(self, name):
        sequencer_model, _ = SequencerModel.objects.get_or_create(model="NextSeq",
                                                                  data_naming_convention=DataGeneration.MISEQ)
        sequencer, _ = Sequencer.objects.get_or_create(name="NS500", sequencer_model=sequencer_model)
        sequencing_run = SequencingRun.objects.create(name=name, sequencer=sequencer)
        sample_sheet = SampleSheet.objects.create(sequencing_run=sequencing_run, hash=name, path=f"/data/{name}.csv")
        SequencingRunCurrentSampleSheet.objects.create(sequencing_run=sequencing_run, sample_sheet=sample_sheet)
        arms = {}
        for i, (extraction, callers) in enumerate([(self.dna, DNA_CALLERS), (self.rna, RNA_CALLERS)], start=1):
            sequencing_sample = SequencingSample.objects.create(sample_sheet=sample_sheet, sample_id=f"SMP{i}",
                                                                sample_name=f"SMP{i}", sample_number=i,
                                                                barcode="ACGT", extraction=extraction)
            arms[extraction] = [SingleSampleVCF.objects.create(path=f"/data/{name}/SMP{i}.{caller}.vcf",
                                                               sequencing_run=sequencing_run,
                                                               sequencing_sample=sequencing_sample,
                                                               variant_caller=self.callers[caller])
                                for caller in callers]
        return sequencing_run, arms

    def _import(self, extraction, single_sample_vcf=None, import_status=ImportStatus.SUCCESS) -> Sample:
        name = single_sample_vcf.path if single_sample_vcf else f"upload_{extraction and extraction.reference_id}"
        vcf = VCF.objects.create(name=name, date=timezone.now(), user=self.user, genotype_samples=1,
                                 genome_build=self.grch37, import_status=import_status, source="DRAGEN")
        sample = Sample.objects.create(vcf=vcf, name=name, import_status=import_status, extraction=extraction)
        assign_permission_to_user_and_groups(self.user, vcf)
        assign_permission_to_user_and_groups(self.user, sample)
        cohort = Cohort.objects.create(name=name, user=self.user, vcf=vcf, genome_build=self.grch37,
                                       import_status=ImportStatus.SUCCESS)
        CohortSample.objects.create(cohort=cohort, sample=sample, cohort_genotype_packed_field_index=0, sort_order=0)
        CohortGenotypeCollection.objects.create(cohort=cohort, cohort_version=cohort.version, num_samples=1)
        if single_sample_vcf:
            file_upload = FileUpload.objects.create(path=name, name=name, user=self.user, import_source="S")
            uploaded_vcf = UploadedVCF.objects.create(file_upload=file_upload, vcf=vcf)
            BackendVCF.objects.create(uploaded_vcf=uploaded_vcf, single_sample_vcf=single_sample_vcf)
            VCFFromSequencingRun.objects.create(vcf=vcf, sequencing_run=single_sample_vcf.sequencing_run,
                                                variant_caller=single_sample_vcf.variant_caller)
        auto_launch_group_analyses_for_sample(self.user, sample)
        return sample

    def _runs(self):
        return GroupAnalysisTemplateRun.objects.filter(specimen=self.specimen)


class AutoLaunchGroupTestCase(AutoLaunchTestBase):

    def test_launches_once_every_registered_file_has_imported(self):
        _, arms = self._make_run("RUN1")
        for single_sample_vcf in arms[self.dna]:
            self._import(self.dna, single_sample_vcf)
        self.assertFalse(self._runs().exists(), "RNA arm still to come")

        self._import(self.rna, arms[self.rna][0])
        run = self._runs().get()
        self.assertEqual(run.analysis_template_run.trigger,
                         AnalysisTemplateRunTrigger.ALL_LINKED_SEQUENCING_SAMPLES_IMPORTED)
        self.assertEqual(run.analysis_template_run.analysis.name, f"{self.template.active} for {self.specimen}")
        node = SampleNode.objects.get(analysis=run.analysis_template_run.analysis)
        self.assertEqual(node.specimen, self.specimen)
        self.assertEqual(len(node.get_sample_group().samples), 3)

    def test_failed_import_stops_the_launch(self):
        _, arms = self._make_run("RUN1")
        self._import(self.dna, arms[self.dna][0])
        self._import(self.dna, arms[self.dna][1], import_status=ImportStatus.ERROR)
        with self.assertLogs(level="WARNING"):  # create_event logs at its severity
            self._import(self.rna, arms[self.rna][0])
        self.assertFalse(self._runs().exists())
        self.assertTrue(Event.objects.filter(name="auto_launch_analysis_template").exists())

    def test_each_sequencing_run_gets_its_own_analysis(self):
        for run_name in ("RUN1", "RUN2"):
            _, arms = self._make_run(run_name)
            for extraction, single_sample_vcfs in arms.items():
                for single_sample_vcf in single_sample_vcfs:
                    self._import(extraction, single_sample_vcf)
        self.assertEqual(self._runs().count(), 2)

    def test_without_sequencing_records_launches_on_first_vcf(self):
        self._import(self.dna)
        self._import(self.rna)
        run = self._runs().get()
        self.assertEqual(run.analysis_template_run.trigger, AnalysisTemplateRunTrigger.VCF_IMPORT)

    def test_sample_linked_after_import_queues_a_retry(self):
        sample = self._import(None)
        self.assertFalse(self._runs().exists())
        with patch("analysis.signals.signal_handlers.auto_run_group_analyses_for_sample") as task:
            with self.captureOnCommitCallbacks(execute=True):
                sample.name = "renamed"
                sample.save()
            task.si.assert_not_called()

            with self.captureOnCommitCallbacks(execute=True):
                sample.extraction = self.dna
                sample.save()
            task.si.assert_called_once_with(sample.pk)


class AutoLaunchMatchingTestCase(AutoLaunchTestBase):

    def test_variant_caller_and_source_regexes(self):
        AutoLaunchAnalysisTemplate.objects.filter(template=self.template).update(variant_caller_regex="TSO_fusion",
                                                                                vcf_source_regex="DRAGEN")
        subjects = {
            (True, True): AutoLaunchSubject(variant_caller_name="TSO_fusion_caller", vcf_source="DRAGEN 4.2"),
            (False, True): AutoLaunchSubject(variant_caller_name="TSO_dragen_caller", vcf_source="DRAGEN 4.2"),
            (False, False): AutoLaunchSubject(variant_caller_name=None, vcf_source=None),
        }
        for (caller_match, source_match), subject in subjects.items():
            match = get_auto_launch_analysis_template_matches(self.user, subject)[0]
            self.assertEqual((match.variant_caller_match, match.vcf_source_match), (caller_match, source_match))

    def test_launch_level_needs_a_template_variable(self):
        auto_launch = AutoLaunchAnalysisTemplate(template=self.template)
        form = AutoLaunchAnalysisTemplateForm({"launch_level": SampleSourceLevel.PATIENT}, instance=auto_launch)
        self.assertIn("launch_level", form.errors)
        form = AutoLaunchAnalysisTemplateForm({"launch_level": SampleSourceLevel.SPECIMEN}, instance=auto_launch)
        self.assertTrue(form.is_valid(), form.errors)

    def test_variant_caller_regex_shows_registered_matches(self):
        auto_launch = AutoLaunchAnalysisTemplate(template=self.template)
        cases = {
            "": None,
            "TSO_dragen": ["TSO_dragen_caller", "TSO_dragen_cnv_caller"],
            "TSO_dargen": [],
        }
        for regex, matches in cases.items():
            form = AutoLaunchAnalysisTemplateForm({"variant_caller_regex": regex}, instance=auto_launch)
            self.assertEqual(form.variant_caller_matches, matches, regex)
