from typing import List

from django.conf import settings
from django.contrib.auth.models import User, Group
from django.contrib.postgres.fields import DecimalRangeField
from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE, PROTECT
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver
from django.urls.base import reverse
from django_dag.models import node_factory, edge_factory
from django_extensions.db.models import TimeStampedModel
import logging
import os
import re
import shutil

from genes.models import CanonicalTranscriptCollection, GeneListCategory, \
    CustomTextGeneList, GeneList, GeneCoverageCollection, Transcript, GeneSymbol, SampleGeneList
from library.enums.log_level import LogLevel
from library.file_utils import name_from_filename, remove_gz_if_exists
from library.log_utils import get_traceback, log_traceback
from library.utils import sorted_nicely
from patients.models import FakeData
from seqauto.illumina import illuminate_report, illumina_sequencers
from seqauto.illumina.illumina_sequencers import SEQUENCING_RUN_REGEX
from seqauto.models_enums import DataGeneration, SequencerRead, PairedEnd, \
    DataState, SequencingFileType, JobScriptStatus, SeqAutoRunStatus, EnrichmentKitType
from seqauto.qc.exec_summary import load_exec_summary
from seqauto.qc.fastqc_parser import read_fastqc_data
from seqauto.qc.flag_stats import load_flagstats
from seqauto.qc.qc_utils import meta_data_file
from snpdb.models import Manufacturer, VCF, Sample, SoftwareVersion, Wiki, LabProject, \
    GenomicIntervalsCollection, GenomeBuild
from snpdb.models.models_enums import ImportStatus, ImportSource
from variantgrid.celery import app


class SeqAutoRun(TimeStampedModel):
    """ Represents a scan and loading of data from backend system """

    status = models.CharField(max_length=1, choices=SeqAutoRunStatus.choices, default=SeqAutoRunStatus.CREATED)
    task_id = models.CharField(max_length=36, null=True)
    scan_start = models.DateTimeField(null=True)
    create_models_start = models.DateTimeField(null=True)
    scripts_and_jobs_start = models.DateTimeField(null=True)
    job_launch_script_filename = models.TextField(null=True)
    finish_date = models.DateTimeField(null=True)
    error_exception = models.TextField(null=True)
    fake_data = models.ForeignKey(FakeData, null=True, on_delete=CASCADE)

    def save(self, **kwargs):
        self.status = self.get_status()
        super().save(**kwargs)

    def get_status(self):
        status = SeqAutoRunStatus.CREATED
        if self.error_exception:
            status = SeqAutoRunStatus.ERROR
        else:
            if self.scan_start:
                status = SeqAutoRunStatus.SCANNING_FILES
            if self.create_models_start:
                status = SeqAutoRunStatus.CREATE_MODELS
            if self.scripts_and_jobs_start:
                status = SeqAutoRunStatus.SCRIPTS_AND_JOBS
            if self.finish_date:
                status = SeqAutoRunStatus.FINISHED
        return status

    def get_scan_resources_dir(self):
        return os.path.join(settings.SEQAUTO_SCAN_RESOURCES_DIR, "seqauto_run_%00d" % self.pk)

    def get_job_scripts_dir(self):
        return os.path.join(settings.SEQAUTO_JOB_SCRIPTS_DIR, "seqauto_run_%00d" % self.pk)

    def remove_scan_resources_dir(self):
        scan_resources_dir = self.get_scan_resources_dir()
        logging.info("*** Deleting files for SeqAutoRun %d - '%s'", self.pk, scan_resources_dir)
        shutil.rmtree(scan_resources_dir)

    @staticmethod
    def get_last_success_datetime():
        successful_runs = SeqAutoRun.objects.filter(status=SeqAutoRunStatus.FINISHED)
        most_recently_finished_seqauto_run = successful_runs.order_by("finish_date").last()
        if most_recently_finished_seqauto_run:
            finish_date = most_recently_finished_seqauto_run.finish_date
        else:
            finish_date = None
        return finish_date

    def __str__(self):
        return f"AnnotationRun: {self.created} ({self.status})"


class SeqAutoRunEvent(models.Model):
    """ These are messages created by SeqAutoRun """
    seqauto_run = models.ForeignKey(SeqAutoRun, on_delete=CASCADE)
    severity = models.CharField(max_length=1, choices=LogLevel.CHOICES)
    file_type = models.CharField(max_length=1, choices=SequencingFileType.choices, null=True)
    message = models.TextField(null=True)

    def __str__(self):
        return f"{self.seqauto_run} {self.get_severity_display()} {self.file_type}: {self.message}"


class SequencerModel(models.Model):
    model = models.TextField(primary_key=True)
    manufacturer = models.ForeignKey(Manufacturer, null=True, on_delete=CASCADE)
    data_naming_convention = models.CharField(max_length=1, choices=DataGeneration.choices)

    @property
    def css_class(self):
        return self.model.replace(' ', '-').lower()

    def __str__(self):
        return f"{self.manufacturer} {self.model}"


class Sequencer(models.Model):
    name = models.TextField(primary_key=True)
    sequencer_model = models.ForeignKey(SequencerModel, on_delete=CASCADE)

    @staticmethod
    def get_sequencers_dict_by_name():
        sequencers = {}
        for sequencer in Sequencer.objects.all():
            sequencers[sequencer.pk] = sequencer
        return sequencers

    @staticmethod
    def get_or_create_sequencer(name, existing_sequencers=None):
        if existing_sequencers is None:
            existing_sequencers = Sequencer.get_sequencers_dict_by_name()

        sequencer = existing_sequencers.get(name)
        if not sequencer:
            model_name = illumina_sequencers.get_sequencer_model_from_name(name)
            if model_name:
                illumina, _ = Manufacturer.objects.get_or_create(name='Illumina')
                try:
                    sequencer_model = SequencerModel.objects.get(model=model_name,
                                                                 manufacturer=illumina)
                except:
                    if 'HiSeq' in model_name:
                        data_naming_convention = DataGeneration.HISEQ
                    else:
                        data_naming_convention = DataGeneration.MISEQ
                    sequencer_model = SequencerModel.objects.create(model=model_name,
                                                                    data_naming_convention=data_naming_convention,
                                                                    manufacturer=illumina)

                sequencer = Sequencer.objects.create(name=name,
                                                     sequencer_model=sequencer_model)
                existing_sequencers[sequencer.name] = sequencer  # Add in case this is being used elsewhere
            else:
                msg = f"Can't recognise Illumina model from name '{name}' and don't know how to handle non-Illumina sequencers yet"
                raise ValueError(msg)

        return sequencer

    def get_absolute_url(self):
        return reverse('view_sequencer', kwargs={'pk': self.pk})

    def __str__(self):
        return f"{self.name} ({self.sequencer_model})"


class EnrichmentKit(models.Model):
    """ A lab method to enrich a sample (eg Capture Panel or Amplicon etc) """
    name = models.TextField()
    version = models.IntegerField(default=1)
    enrichment_kit_type = models.CharField(max_length=1, choices=EnrichmentKitType.choices, null=True)
    manufacturer = models.ForeignKey(Manufacturer, null=True, on_delete=CASCADE)
    bed_file = models.TextField(null=True, blank=True)  # Original manifest bed
    genomic_intervals = models.ForeignKey(GenomicIntervalsCollection, null=True, blank=True, on_delete=SET_NULL)
    gene_list = models.ForeignKey(GeneList, null=True, on_delete=PROTECT)
    canonical_transcript_collection = models.ForeignKey(CanonicalTranscriptCollection, null=True, on_delete=SET_NULL)
    obsolete = models.BooleanField(default=False)

    @staticmethod
    def get_latest_version(name):
        return EnrichmentKit.objects.filter(name=name).order_by("version").last()

    @staticmethod
    def get_enrichment_kits(enrichment_kit_kwargs_list):
        enrichment_kits = []
        for enrichment_kit_kwargs in enrichment_kit_kwargs_list:
            # There can be more than 1 per enrichment_kit name
            for enrichment_kit in EnrichmentKit.objects.filter(**enrichment_kit_kwargs):
                enrichment_kits.append(enrichment_kit)
        return enrichment_kits

    def get_gold_sequencing_runs_qs(self):
        return self.sequencingrun_set.filter(gold_standard=True)

    def get_absolute_url(self):
        return reverse('view_enrichment_kit', kwargs={"pk": self.pk})

    def __str__(self):
        name = self.name
        if self.version > 1:
            name += f" (version {self.version})"
        return name


class Library(models.Model):
    name = models.TextField(primary_key=True)
    description = models.TextField(null=True, blank=True)
    manufacturer = models.ForeignKey(Manufacturer, null=True, on_delete=CASCADE)

    def get_absolute_url(self):
        return reverse('view_library', kwargs={'pk': self.pk})

    def __str__(self):
        return ' '.join(map(str, filter(lambda x: x, [self.name, self.description, self.manufacturer])))


class Assay(models.Model):
    library = models.ForeignKey(Library, on_delete=CASCADE)
    sequencer = models.ForeignKey(Sequencer, on_delete=CASCADE)
    enrichment_kit = models.ForeignKey(EnrichmentKit, on_delete=CASCADE)  # Genes captured

    def get_absolute_url(self):
        return reverse('view_assay', kwargs={'pk': self.pk})

    def __str__(self):
        return f"Seq: {self.sequencer}, Lib: {self.library}, EnrichmentKit: {self.enrichment_kit}"


class VariantCaller(SoftwareVersion):
    run_params = models.TextField(null=True, blank=True)

    def __str__(self):
        description = f"{self.name} (v. {self.version})"
        if self.run_params:
            description += f": {self.run_params}"
        return description

    def get_absolute_url(self):
        return reverse('view_variant_caller', kwargs={'pk': self.pk})


class Aligner(SoftwareVersion):

    def __str__(self):
        return f"{self.name} (v. {self.version})"

    def get_absolute_url(self):
        return reverse('view_aligner', kwargs={'pk': self.pk})


class VariantCallingPipeline(SoftwareVersion):
    description = models.TextField(null=True, blank=True)
    aligner = models.ForeignKey(Aligner, on_delete=CASCADE)
    variant_caller = models.ForeignKey(VariantCaller, on_delete=CASCADE)
    other_details = models.TextField(null=True, blank=True)

    def get_absolute_url(self):
        return reverse('view_variant_calling_pipeline', kwargs={'pk': self.pk})

    def __str__(self):
        return f"{self.description} (Algn: {self.aligner}, V. Caller: {self.variant_caller})"

# Software Pipeline Models replace VariantCallingPipeline and use the acyclic
# graph module from django_dag


class SoftwarePipeline(SoftwareVersion):
    description = models.TextField(null=True, blank=True)


class SoftwarePipelineNode(node_factory('SoftwarePipelineEdge')):
    softwarepipeline = models.ForeignKey(SoftwarePipeline, on_delete=CASCADE)
    name = models.TextField()
    version = models.TextField()
    parameters = models.TextField()
    description = models.TextField()


class SoftwarePipelineEdge(edge_factory(SoftwarePipelineNode, concrete=False)):
    name = models.CharField(max_length=32, blank=True, null=True)


# ExperimentManager manages the Experiment objects
# to ensure that experiment names are cleaned before creation or updating
class ExperimentManager(models.Manager):

    @staticmethod
    def fix_kwargs(kwargs):
        try:
            old_name = kwargs["name"]
            name = Experiment.clean_experiment_name(old_name)
            kwargs["name"] = name
        except:
            pass
        return kwargs

    def get(self, *args, **kwargs):
        self.fix_kwargs(kwargs)
        return super().get(*args, **kwargs)

    def get_or_create(self, default=None, **kwargs):
        self.fix_kwargs(kwargs)
        return super().get_or_create(default, **kwargs)


class Experiment(models.Model):
    name = models.TextField(primary_key=True)
    created = models.DateTimeField(auto_now_add=True)
    objects = ExperimentManager()

    @staticmethod
    def clean_experiment_name(experiment_name):
        experiment_name = experiment_name.upper()
        experiment_name = experiment_name.replace("-", "_")
        experiment_name = experiment_name.replace(" ", "_")

        # Remove RPT off the end...
        experiment_name = re.sub("_RPT$", "", experiment_name)
        experiment_name = re.sub("RPT$", "", experiment_name)
        return experiment_name

    def save(self, **kwargs):
        old_name = self.name
        self.name = Experiment.clean_experiment_name(old_name)
        return super().save(**kwargs)

    def can_write(self, user):
        """ can't delete once you've linked to SequencingRun """
        return user.is_superuser and not self.sequencingrun_set.exists()

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('view_experiment', kwargs={"experiment_id": self.pk})


class SequencingRun(models.Model):
    """ Represents a flowcell (or other technology with multiple sequencing samples) """
    name = models.TextField(primary_key=True)
    sequencer = models.ForeignKey(Sequencer, on_delete=CASCADE)
    path = models.TextField()
    gold_standard = models.BooleanField(default=False)
    bad = models.BooleanField(default=False)
    hidden = models.BooleanField(default=False)
    ready = models.BooleanField(default=False)
    experiment = models.ForeignKey(Experiment, null=True, on_delete=SET_NULL)
    enrichment_kit = models.ForeignKey(EnrichmentKit, null=True, on_delete=CASCADE)  # Sequencing Run can be all one enrichment_kit, or SequencingSample can have own enrichment_kits
    data_state = models.CharField(max_length=1, choices=DataState.choices)  # This is for the SequencingRun directory
    has_basecalls = models.BooleanField(default=False)
    has_interop = models.BooleanField(default=False)  # Quality, Index and Tile
    fake_data = models.ForeignKey(FakeData, null=True, on_delete=CASCADE)

    def save(self, force_insert=False, force_update=False, using=None,
             update_fields=None):
        self.ready = not bool(self.get_errors())
        super().save()

    def get_current_sample_sheet(self):
        try:
            return self.sequencingruncurrentsamplesheet.sample_sheet
        except:
            logging.info("Can't find current sample sheet for %s", self.path)
            raise

    def get_params(self):
        """ This allows chaining down names etc - in case a level changes it, will cascade down """
        params = {"sequencing_run": self.name,
                  "sequencing_run_dir": self.path}

        # TAU rename the sequencing run dir with enrichment kit at the end - need to clean it
        original_sequencing_run = self.name
        if m := re.match(SEQUENCING_RUN_REGEX, self.name):
            original_sequencing_run = m.group(0)
        params["original_sequencing_run"] = original_sequencing_run
        if self.enrichment_kit:
            params["enrichment_kit"] = self.enrichment_kit.name

        if self.experiment:
            params["experiment"] = self.experiment.name

            # TODO: TAU hack - not sure how we're supposed to get this
            TAU_LAB_IDS = {"idt_exome": "Exome", "idt_gmp_focus": "GMPFocus"}
            if self.enrichment_kit:
                if nice_name := TAU_LAB_IDS.get(self.enrichment_kit.name):
                    experiment_parts = self.experiment.name.split("_")
                    if len(experiment_parts) >= 3:  # Could be GMPFOCUS_21_011 or GMPFOCUS_21_011_FFPE
                        lab_run_id = f"{nice_name}_{experiment_parts[1]}_{experiment_parts[2]}"
                        params["lab_run_id"] = lab_run_id
        return params

    def get_enrichment_kit_name(self):
        if self.enrichment_kit:
            enrichment_kit_name = self.enrichment_kit.name
        else:
            enrichment_kit_name = 'Unknown EnrichmentKit'
        return enrichment_kit_name

    def check_basecalls_dir(self):
        basecalls_dir = os.path.join(self.path, "Data", "Intensities", "BaseCalls")
        has_basecall_data = False
        if os.path.exists(basecalls_dir):
            for f in os.listdir(basecalls_dir):
                if os.path.isdir(f) and f.startswith("L00"):
                    has_basecall_data = True
        return has_basecall_data

    def check_interop_dir(self):
        interop_dir = os.path.join(self.path, "InterOp")
        has_quality = os.path.exists(os.path.join(interop_dir, "QMetricsOut.bin"))
        has_index = os.path.exists(os.path.join(interop_dir, "IndexMetricsOut.bin"))
        has_tile = os.path.exists(os.path.join(interop_dir, "TileMetricsOut.bin"))
        return has_quality and has_index and has_tile

    def can_basecall(self):
        return self.has_basecalls

    def can_generate_qc(self):
        return self.data_state != DataState.DELETED and self.has_interop

    def get_old_sample_sheets(self):
        current_sample_sheet = self.get_current_sample_sheet()
        return SampleSheet.objects.filter(sequencing_run=self).exclude(pk=current_sample_sheet.pk)

    def get_errors(self) -> List[str]:
        errors = []
        if self.is_data_out_of_date_from_current_sample_sheet:
            errors.append(f"SampleSheet has changed, please confirm in 'Admin' tab")
        try:
            illumina_qc = self.sequencingruncurrentsamplesheet.sample_sheet.illuminaflowcellqc
            if illumina_qc.data_state != DataState.COMPLETE:
                errors.append(f"QC data state={illumina_qc.get_data_state_display()}")
        except SequencingRunCurrentSampleSheet.DoesNotExist:
            errors.append("No Current SampleSheet set")
        except:
            errors.append("QC not loaded")
        return errors

    def get_warnings(self) -> List[str]:
        warnings = []
        return warnings

    @property
    def is_data_out_of_date_from_current_sample_sheet(self):
        current_sample_sheet = self.get_current_sample_sheet()
        try:
            combo = current_sample_sheet.samplesheetcombinedvcffile_set.get()
            if combo.needs_to_be_linked():
                return True
        except:
            pass

        old_sample_sheets = self.get_old_sample_sheets()
        illuminate_qc = IlluminaFlowcellQC.objects.filter(sample_sheet=current_sample_sheet)
        old_illuminate_qc = IlluminaFlowcellQC.objects.filter(sample_sheet__in=old_sample_sheets)
        old_sample_sheet_combined_vcf = SampleSheetCombinedVCFFile.objects.filter(sample_sheet__in=old_sample_sheets)
        old_sample_links = SampleFromSequencingSample.objects.filter(sequencing_sample__sample_sheet__in=old_sample_sheets)
        old_unaligned_reads = UnalignedReads.get_for_old_sample_sheets(self)
        return any([not illuminate_qc.exists() and old_illuminate_qc.exists(),
                    old_sample_sheet_combined_vcf.exists(),
                    old_unaligned_reads.exists(),
                    old_sample_links.exists()])

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        return reverse('view_sequencing_run', kwargs={"sequencing_run_id": self.pk})


class SequencingRunWiki(Wiki):
    sequencing_run = models.OneToOneField(SequencingRun, on_delete=CASCADE)


class SequencingRunModification(models.Model):
    """ To be able to associate a sequencing_run with an event (or just general message) """

    created = models.DateTimeField(auto_now_add=True)
    sequencing_run = models.ForeignKey(SequencingRun, on_delete=CASCADE)
    user = models.ForeignKey(User, null=True, on_delete=SET_NULL)
    message = models.TextField(null=True)  # Filled in by users
    seqauto_run_event = models.ForeignKey(SeqAutoRunEvent, null=True, on_delete=CASCADE)  # Set if done as an event


class SeqAutoFile(models.Model):
    path = models.TextField()
    # last modified - Stores stat st_mtime - time of last modification - only used for classes that can reload
    file_last_modified = models.FloatField(default=0.0)
    data_state = models.CharField(max_length=1, choices=DataState.choices)
    error_exception = models.TextField(null=True)
    #hash = models.TextField() ???

    class Meta:
        abstract = True


class SampleSheetPropertiesMixin:
    """ Needs a "sample_sheet" property """

    @property
    def sequencing_run(self):
        return self.sample_sheet.sequencing_run


class SequencingSamplePropertiesMixin(SampleSheetPropertiesMixin):
    """ Needs a "sequencing_sample" property """

    @property
    def sample_sheet(self):
        return self.sequencing_sample.sample_sheet


class SampleSheet(models.Model):
    SAMPLE_SHEET = "SampleSheet.csv"
    # We should assume that these will be changed, and so need to know...
    hash = models.TextField()
    sequencing_run = models.ForeignKey(SequencingRun, on_delete=CASCADE)
    path = models.TextField()
    date = models.DateField(null=True)

    def get_params(self):
        params = self.sequencing_run.get_params()
        params["sample_sheet"] = self.path
        return params

    def get_sequencing_samples_by_name(self):
        sequencing_samples_by_name = {}
        for ss in self.sequencingsample_set.all():
            sequencing_samples_by_name[ss.sample_name] = ss
        return sequencing_samples_by_name

    def get_sorted_sequencing_samples(self):
        samples_by_name = self.get_sequencing_samples_by_name()
        return [samples_by_name[name] for name in sorted_nicely(samples_by_name)]

    def get_sample_enrichment_kits(self):
        qs = self.sequencingsample_set.filter(enrichment_kit__isnull=False)
        return EnrichmentKit.objects.filter(pk__in=qs.values_list("enrichment_kit", flat=True))

    @staticmethod
    def get_path_from_sequencing_run(sequencing_run):
        sample_sheet = os.path.join(sequencing_run.path, SampleSheet.SAMPLE_SHEET)
        return os.path.abspath(sample_sheet)

    def get_version_string(self):
        return f"{self.hash}/{self.date}"

    def __str__(self):
        return self.path


class SequencingRunCurrentSampleSheet(models.Model):
    """ This is a way to "walk" via relations from sequencing run to the latest sample sheet """
    sequencing_run = models.OneToOneField(SequencingRun, on_delete=CASCADE)
    sample_sheet = models.OneToOneField(SampleSheet, on_delete=CASCADE)

    def __str__(self):
        return f"{self.sequencing_run}, {self.sample_sheet}"


class SequencingSample(models.Model):
    sample_sheet = models.ForeignKey(SampleSheet, on_delete=CASCADE)
    sample_id = models.TextField()
    # sample_name is used to name files. In MiSeq/NextSeq samplesheet you can add names.
    # For Hiseq and if left empty on MiSeq this will be sample_id
    sample_name = models.TextField(null=True)
    sample_project = models.TextField(null=True)
    sample_number = models.IntegerField()
    lane = models.IntegerField(null=True)
    barcode = models.TextField()
    enrichment_kit = models.ForeignKey(EnrichmentKit, null=True, on_delete=CASCADE)
    is_control = models.BooleanField(default=False)
    failed = models.BooleanField(default=False)
    automatically_process = models.BooleanField(default=True)

    def get_params(self):
        params = self.sample_sheet.get_params()
        params.update({"lane": self.lane or 1,
                       "sample_number": self.sample_number,
                       "sample_id": self.sample_id,
                       "sample_name": self.sample_name,
                       "sample_project": self.sample_project or '',
                       "sample_name_underscores": self.sample_name.replace("-", "_"),
                       "barcode": self.barcode})

        if self.enrichment_kit is not None:
            params['enrichment_kit'] = self.enrichment_kit.name

        sequencer_model = self.sample_sheet.sequencing_run.sequencer.sequencer_model
        if sequencer_model.data_naming_convention == DataGeneration.HISEQ:
            full_sample_name = f"{self.sample_id}_{self.barcode}"
        elif sequencer_model.data_naming_convention == DataGeneration.MISEQ:
            full_sample_name = f"{self.sample_name}_S{self.sample_number}"
        else:
            msg = f"Unknown sequencer_model.data_naming_convention '{sequencer_model.data_naming_convention}'"
            raise ValueError(msg)
        params['full_sample_name'] = full_sample_name
        params['full_sample_name_underscores'] = full_sample_name.replace("-", "_")
        return params

    def get_single_bam(self):
        """ relies on there being only one """
        try:
            unaligned_reads = self.unalignedreads_set.get()
            try:
                bam_file = unaligned_reads.bamfile_set.get()
                return bam_file
            except:
                logging.error("Wasn't exactly 1 bam_file for unaligned_reads %s", unaligned_reads)
        except:
            logging.error("Wasn't exactly 1 unaligned reads for sequencing sample %s", self)
        return None

    def get_single_qc(self):
        bam_file = self.get_single_bam()
        return bam_file.qc_set.get()

    @staticmethod
    def get_current():
        """ Return SequencingSamples that have not been replaced by newer SampleSheet """
        return SequencingSample.objects.filter(sample_sheet__sequencingruncurrentsamplesheet__isnull=False)

    def __str__(self):
        return self.sample_id


class SequencingSampleData(models.Model):
    """ key/values set from settings.SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS """
    sequencing_sample = models.ForeignKey(SequencingSample, on_delete=CASCADE)
    column = models.TextField()
    value = models.TextField(null=True)


# Store the final relations, so it's easier to move between these 2 models
class SampleFromSequencingSample(models.Model):
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    sequencing_sample = models.OneToOneField(SequencingSample, on_delete=CASCADE)

    @property
    def sequencing_run(self):
        return self.sequencing_sample.sample_sheet.sequencing_run


class VCFFromSequencingRun(models.Model):
    vcf = models.OneToOneField(VCF, on_delete=CASCADE)
    sequencing_run = models.OneToOneField(SequencingRun, on_delete=CASCADE)


class IlluminaFlowcellQC(SeqAutoFile, SampleSheetPropertiesMixin):
    sample_sheet = models.OneToOneField(SampleSheet, on_delete=CASCADE)
    mean_cluster_density = models.IntegerField(null=True)
    mean_pf_cluster_density = models.IntegerField(null=True)
    total_clusters = models.IntegerField(null=True)
    total_pf_clusters = models.IntegerField(null=True)
    percentage_of_clusters_pf = models.FloatField(null=True)
    aligned_to_phix = models.FloatField(null=True)

    def get_params(self):
        params = self.sample_sheet.get_params()
        params['illuminate_qc'] = self.path
        return params

    def load_from_file(self, seqauto_run, **kwargs):
        illuminate_report.load_from_file(seqauto_run, self)

    @staticmethod
    def get_path_from_sequencing_run(sequencing_run):
        illuminate_path = os.path.join(sequencing_run.path, "Illuminate", 'illuminate_report.txt')
        return os.path.abspath(illuminate_path)

    @staticmethod
    def get_sequencing_run_path():
        return 'sample_sheet__sequencing_run'

    def __str__(self):
        return self.path


class ReadQ30(models.Model):
    illumina_flowcell_qc = models.ForeignKey(IlluminaFlowcellQC, on_delete=CASCADE)
    sequencer_read_id = models.IntegerField()  # Eg HiSeq= [R1,Index,R2], NextSeq/MiSeq=[R1,Index1,Index2,R2]
    read = models.CharField(max_length=2, choices=SequencerRead.choices)
    percent = models.FloatField()
    is_index = models.BooleanField(default=False)

    def __str__(self):
        index = " (Index)" if self.is_index else ''
        return f"Read: {self.sequencer_read_id}{index}. {self.read} {self.percent:.2f}%"


class IlluminaIndexQC(models.Model):
    illumina_flowcell_qc = models.ForeignKey(IlluminaFlowcellQC, on_delete=CASCADE)
    index = models.TextField()
    project = models.TextField()
    name = models.TextField()
    reads = models.IntegerField()


class Fastq(SeqAutoFile, SequencingSamplePropertiesMixin):
    sequencing_sample = models.ForeignKey(SequencingSample, on_delete=CASCADE)
    name = models.TextField()  # from path
    read = models.CharField(max_length=2, choices=PairedEnd.choices)

    def get_common_pair_name_and_read(self):
        name = os.path.basename(self.path)
        regex_pattern = "(.+)_(R[12])(_001)?.fastq.gz"

        m = re.match(regex_pattern, name)
        if not m:
            msg = f"FastQ path '{self.path}' does not match regex pattern '{regex_pattern}'"
            raise ValueError(msg)
        pair_name = m.group(1)
        read = m.group(2)
        return pair_name, read

    def get_params(self):
        params = self.sequencing_sample.get_params()
        params['fastq'] = self.path
        return params

    @staticmethod
    def get_pair_paths_from_sequencing_sample(sequencing_sample):
        """ The 1st value returned must be what will be generated given current pipelines (ie old naming conventions come after) """
        sequencer_model = sequencing_sample.sample_sheet.sequencing_run.sequencer.sequencer_model
        patterns = ["%(full_sample_name)s_L%(lane)03d_R%(read)d_001.fastq.gz"]

        if sequencer_model.data_naming_convention == DataGeneration.HISEQ:
            patterns.append("%(sample_id)s_R%(read)d.fastq.gz")  # SACGF style naming convention
        elif sequencer_model.data_naming_convention == DataGeneration.MISEQ:
            patterns.extend(["%(sample_name)s_S%(sample_number)d_R%(read)d_001.fastq.gz",
                             "%(sample_name_underscores)s_S%(sample_number)d_L%(lane)03d_R%(read)d_001.fastq.gz",
                             "%(sample_name_underscores)s_S%(sample_number)d_R%(read)d_001.fastq.gz"])
        else:
            msg = f"Unknown sequencer_model.data_naming_convention '{sequencer_model.data_naming_convention}'"
            raise ValueError(msg)

        params = sequencing_sample.get_params()
        r1_params = {"read": 1}
        r1_params.update(params)
        r2_params = {"read": 2}
        r2_params.update(params)

        unaligned_dir_patterns = []
        if sequencing_sample.sample_project:
            vcf_subdirs = ["%(sequencing_run_dir)s/Data/Intensities/BaseCalls/%(sample_project)s",
                           "%(sequencing_run_dir)s/Data/Intensities/BaseCalls/%(sample_project)s/%(sample_id)s"]
            unaligned_dir_patterns.extend(vcf_subdirs)
        unaligned_dir_patterns.append("%(sequencing_run_dir)s/Data/Intensities/BaseCalls")

        pair_paths = []
        for unaligned_dir_pattern in unaligned_dir_patterns:
            fastq_dir = unaligned_dir_pattern % params

            for pattern in patterns:
                name1 = pattern % r1_params
                name2 = pattern % r2_params

                fastq_r1 = os.path.abspath(os.path.join(fastq_dir, name1))
                fastq_r2 = os.path.abspath(os.path.join(fastq_dir, name2))
                pair_paths.append((fastq_r1, fastq_r2))

        return pair_paths

    def __str__(self):
        return f"FastQ {self.name} ({self.read}) from sample {self.sequencing_sample}"


class FastQC(SeqAutoFile, SequencingSamplePropertiesMixin):
    fastq = models.OneToOneField(Fastq, on_delete=CASCADE)
    total_sequences = models.IntegerField(null=True)
    filtered_sequences = models.IntegerField(null=True)
    gc = models.IntegerField(null=True)

    @property
    def sequencing_sample(self):
        return self.fastq.sequencing_sample

    def get_params(self):
        params = self.fastq.get_params()
        return params

    def load_from_file(self, seqauto_run, **kwargs):
        fastqc_data = read_fastqc_data(self.path)
        basic_stats = fastqc_data['Basic Statistics']

        self.total_sequences = basic_stats["Total Sequences"]
        self.filtered_sequences = basic_stats['Sequences flagged as poor quality']
        self.gc = basic_stats["GC"]

    @staticmethod
    def get_path_from_fastq(fastq):
        base_dir = os.path.dirname(fastq.path)
        name = remove_gz_if_exists(fastq.path)
        fastqc_dir = "%s_fastqc" % name_from_filename(name)
        path = os.path.join(base_dir, "FastQC", fastqc_dir, "fastqc_data.txt")
        return os.path.abspath(path)

    def __str__(self):
        return f"FastQC ({self.get_data_state_display()}) for {self.fastq}"


class UnalignedReads(models.Model, SequencingSamplePropertiesMixin):
    """ A way of keeping fastqs together """
    sequencing_sample = models.ForeignKey(SequencingSample, on_delete=CASCADE)
    fastq_r1 = models.ForeignKey(Fastq, related_name='fastq_r1', on_delete=CASCADE)
    fastq_r2 = models.ForeignKey(Fastq, null=True, related_name='fastq_r2', on_delete=CASCADE)

    def get_params(self):
        fastq_params = self.fastq_r1.sequencing_sample.get_params()
        data_naming_convention = self.sequencing_run.sequencer.sequencer_model.data_naming_convention
        if data_naming_convention == DataGeneration.HISEQ:
            aligned_pattern = settings.SEQAUTO_HISEQ_ALIGNED_PATTERN
        elif data_naming_convention == DataGeneration.MISEQ:
            aligned_pattern = settings.SEQAUTO_MISEQ_ALIGNED_PATTERN
        else:
            msg = f"Unknown data_naming_convertion: {data_naming_convention}"
            raise ValueError(msg)
        params = {"aligned_pattern": aligned_pattern % fastq_params,
                  "read_1": self.fastq_r1.path}
        if self.fastq_r2:
            params["read_2"] = self.fastq_r2.path

        params.update(fastq_params)
        return params

    @property
    def data_state(self):
        r1_ds = self.fastq_r1.data_state

        if self.fastq_r2 is None:
            return r1_ds

        r2_ds = self.fastq_r2.data_state

        if r1_ds == r2_ds:  # Same - easy
            return r1_ds
        OVERWRITE_OTHER = [DataState.ERROR, DataState.SKIPPED, DataState.DELETED, DataState.NON_EXISTENT]
        for o in OVERWRITE_OTHER:
            if o in (r1_ds, r2_ds):
                return o

        msg = f"Unaligned reads {self} don't know how to get data_state from r1/r2"
        raise ValueError(msg)

    @staticmethod
    def get_for_old_sample_sheets(sequencing_run):
        old_sample_sheets = sequencing_run.get_old_sample_sheets()
        return UnalignedReads.objects.filter(sequencing_sample__sample_sheet__in=old_sample_sheets)

    def __str__(self):
        if self.fastq_r1.data_state == self.fastq_r2.data_state:
            data_state = self.fastq_r1.get_data_state_display()
        else:
            data_state = f"R1: {self.fastq_r1.get_data_state_display()}, R2: {self.fastq_r2.get_data_state_display()}"
        return f"UnalignedReads from sample {self.sequencing_sample} ({data_state})"


class BamFile(SeqAutoFile, SequencingSamplePropertiesMixin):
    unaligned_reads = models.ForeignKey(UnalignedReads, null=True, on_delete=CASCADE)
    name = models.TextField()
    aligner = models.ForeignKey(Aligner, on_delete=CASCADE)

    def get_params(self):
        params = self.unaligned_reads.get_params()
        params['bam'] = self.get_path_from_unaligned_reads(self.unaligned_reads)
        return params

    @property
    def sequencing_sample(self):
        return self.unaligned_reads.sequencing_sample

    @staticmethod
    def get_path_from_unaligned_reads(unaligned_reads):
        params = unaligned_reads.get_params()
        pattern = os.path.join(settings.SEQAUTO_ALIGNED_DIR_PATTERN, settings.SEQAUTO_BAM_PATTERN)
        try:
            filename = pattern % params
        except KeyError as ke:
            if "enrichment_kit" in ke.args:
                logging.error("'enrichment_kit' not set, this is usually done via signal handlers in a custom app")
            logging.error(f"{unaligned_reads} missing: {', '.join(ke.args)}. params={params}")
            raise
        return os.path.abspath(filename)

    @staticmethod
    def get_aligner_from_bam_file(bam_path):
        # TODO: Do properly
        aligner, _ = Aligner.objects.get_or_create(name='fake_aligner', version='0.666')
        return aligner

    def __str__(self):
        return f"BAM: {self.unaligned_reads.sequencing_sample} ({self.get_data_state_display()})"


class Flagstats(SeqAutoFile):
    FLAGSTATS_EXTENSION = ".flagstat.txt"

    bam_file = models.OneToOneField(BamFile, on_delete=CASCADE)
    total = models.IntegerField(null=True)
    read1 = models.IntegerField(null=True)
    read2 = models.IntegerField(null=True)
    mapped = models.IntegerField(null=True)
    properly_paired = models.IntegerField(null=True)

    @property
    def mapped_percent(self):
        return 100.0 * self.mapped / self.total

    @property
    def properly_paired_percent(self):
        return 100.0 * self.properly_paired / self.total

    def get_params(self):
        return self.bam_file.get_params()

    def load_from_file(self, seqauto_run, **kwargs):
        flagstats_data = load_flagstats(self.path)
        self.total = flagstats_data["total"]
        self.read1 = flagstats_data["read1"]
        self.read2 = flagstats_data["read2"]
        self.mapped = flagstats_data["mapped"]
        self.properly_paired = flagstats_data["properly paired"]

    @staticmethod
    def get_path_from_bam_file(bam_file):
        return meta_data_file(bam_file.path, "flagstats/%%s%s" % Flagstats.FLAGSTATS_EXTENSION)

    def __str__(self):
        return f"Flagstats ({self.get_data_state_display()}) for {self.bam_file}"


def get_fake_variant_caller():
    variant_caller, _ = VariantCaller.objects.get_or_create(name='fake_variant_caller',
                                                            version='0.666')
    return variant_caller


def get_seqauto_user():
    user, created = User.objects.get_or_create(username=settings.SEQAUTO_USER)
    if created:
        if settings.SEQAUTO_GROUP:
            seqauto_group = Group.objects.get_or_create(name=settings.SEQAUTO_GROUP)
            user.groups.add(seqauto_group)
    return user


def import_filesystem_vcf(path, name, import_source):
    seqauto_user = get_seqauto_user()

    CELERY_PROCESS_VCF_TASK = 'upload.tasks.vcf.import_vcf_tasks.process_vcf_file_task'
    logging.info("Executing: %s", CELERY_PROCESS_VCF_TASK)
    app.send_task(CELERY_PROCESS_VCF_TASK, args=[path, name, seqauto_user.pk, import_source])  # @UndefinedVariable


class DontAutoLoadException(Exception):
    pass


class VCFFile(SeqAutoFile, SequencingSamplePropertiesMixin):
    """ VCFs from the file system """
    bam_file = models.ForeignKey(BamFile, on_delete=CASCADE)
    variant_caller = models.ForeignKey(VariantCaller, on_delete=CASCADE)

    def get_params(self):
        return self.bam_file.get_params()

    @property
    def name(self):
        return name_from_filename(self.path)

    @property
    def sequencing_sample(self):
        return self.bam_file.sequencing_sample

    @property
    def upload_pipeline(self):
        try:
            up = self.backendvcf.uploaded_vcf.upload_pipeline
        except:
            up = None
        return up

    def load_from_file(self, seqauto_run, **kwargs):
        if not settings.SEQAUTO_IMPORT_VCF:
            raise DontAutoLoadException()

        import_filesystem_vcf(self.path, self.name, ImportSource.SEQAUTO)

    @staticmethod
    def get_path_from_bam(bam):
        pattern = os.path.join(settings.SEQAUTO_VCF_DIR_PATTERN, settings.SEQAUTO_VCF_PATTERN)
        vcf_path = pattern % bam.get_params()
        return os.path.abspath(vcf_path)

    @staticmethod
    def get_variant_caller_from_vcf_file(vcf_path):
        return get_fake_variant_caller()

    def __str__(self):
        return f"VCF {self.name} ({self.get_data_state_display()})"


class SampleSheetCombinedVCFFile(SeqAutoFile, SampleSheetPropertiesMixin):
    sample_sheet = models.ForeignKey(SampleSheet, on_delete=CASCADE)
    variant_caller = models.ForeignKey(VariantCaller, on_delete=CASCADE)

    def get_params(self):
        return self.sample_sheet.get_params()

    @property
    def name(self):
        dir_name = os.path.basename(os.path.dirname(self.path))
        return f"{dir_name} combined"

    @property
    def upload_pipeline(self):
        try:
            up = self.backendvcf.uploaded_vcf.upload_pipeline
        except:
            up = None
        return up

    def load_from_file(self, seqauto_run, **kwargs):
        if not settings.SEQAUTO_IMPORT_COMBO_VCF:
            raise DontAutoLoadException()

        try:
            vcf = VCF.objects.filter(uploadedvcf__uploaded_file__path=self.path)
            if vcf.exists():
                try:
                    from upload.models import UploadedVCF
                    from upload.vcf.vcf_import import create_backend_vcf_links

                    uploaded_vcf = UploadedVCF.objects.filter(uploaded_file__path=self.path).get()
                    create_backend_vcf_links(uploaded_vcf)
                except:
                    log_traceback()

                # Re-linking SequencingRun / backend VCF will be done manually in view_sequencing_run
                # "Assign data to most recent sample sheet" button
                logging.info("Skipping loading of %s as vcf named already exists!", self.path)
                raise DontAutoLoadException()
        except VCF.DoesNotExist:
            pass

        import_filesystem_vcf(self.path, self.name, ImportSource.SEQAUTO)

    def needs_to_be_linked(self):
        try:
            self.backendvcf  # if no exception we can link it
            try:
                self.sample_sheet.sequencing_run.vcffromsequencingrun
            except VCFFromSequencingRun.DoesNotExist:
                return True
        except:
            pass
        return False

    @staticmethod
    def get_path_from_sample_sheet(sample_sheet):
        try:
            pattern = os.path.join(settings.SEQAUTO_VCF_DIR_PATTERN, settings.SEQAUTO_COMBINED_VCF_PATTERN)
            params = sample_sheet.get_params()
            combined_vcf_path = pattern % params
            return os.path.abspath(combined_vcf_path)
        except Exception as e:
            raise ValueError(f"sample_sheet: {sample_sheet}") from e

    @staticmethod
    def get_variant_caller_from_vcf_file(vcf_path):
        return get_fake_variant_caller()

    def __str__(self):
        num_samples = self.sample_sheet.sequencingsample_set.all().count()
        return f"ComboVCF ({self.pk}) for {self.sample_sheet.sequencing_run} ({num_samples} samples)"


class QC(SeqAutoFile, SequencingSamplePropertiesMixin):
    """ This holds together all of the other QC objects """
    bam_file = models.ForeignKey(BamFile, on_delete=CASCADE)
    vcf_file = models.ForeignKey(VCFFile, on_delete=CASCADE)

    @property
    def genome_build(self):
        return GenomeBuild.legacy_build()

    def get_params(self):
        return self.vcf_file.get_params()

    @property
    def qc_dir_path(self):
        return os.path.dirname(self.path)

    @property
    def sequencing_sample(self):
        return self.bam_file.sequencing_sample

    @staticmethod
    def get_path_from_vcf(vcf):
        pattern = os.path.join(settings.SEQAUTO_QC_DIR_PATTERN, settings.SEQAUTO_QC_EXEC_SUMMARY_PATTERN)
        qc_path = pattern % vcf.get_params()
        return os.path.abspath(qc_path)

    @staticmethod
    def get_tsv_path_from_vcf(vcf):
        pattern = os.path.join(settings.SEQAUTO_QC_DIR_PATTERN, settings.SEQAUTO_QC_EXEC_SUMMARY_TSV_PATTERN)
        qc_path = pattern % vcf.get_params()
        return os.path.abspath(qc_path)

    def load_from_file(self, seqauto_run, **kwargs):
        QCExecSummary.load_for_qc(seqauto_run, self, **kwargs)

    def __str__(self):
        return f"QC {name_from_filename(self.path)} ({self.get_data_state_display()})"


class QCGeneList(SeqAutoFile, SequencingSamplePropertiesMixin):
    """ This represents a text file containing genes which will be used for initial pass and QC filters """
    qc = models.ForeignKey(QC, on_delete=CASCADE)
    custom_text_gene_list = models.OneToOneField(CustomTextGeneList, null=True, on_delete=SET_NULL)
    sample_gene_list = models.ForeignKey(SampleGeneList, null=True, on_delete=SET_NULL)

    @property
    def sequencing_sample(self):
        return self.qc.sequencing_sample

    @staticmethod
    def get_path_from_qc(qc):
        pattern = os.path.join(qc.qc_dir_path, settings.SEQAUTO_QC_GENE_LIST_PATTERN)
        return pattern % qc.get_params()

    def load_from_file(self, seqauto_run, **kwargs):
        from genes.custom_text_gene_list import create_custom_text_gene_list
        with open(self.path) as f:
            custom_gene_list_text = f.read()
            custom_text_gene_list = CustomTextGeneList(name=f"QC GeneList for {self.sequencing_sample.sample_name}",
                                                       text=custom_gene_list_text)
            custom_text_gene_list.save()
            seqauto_user = get_seqauto_user()
            create_custom_text_gene_list(custom_text_gene_list, seqauto_user, GeneListCategory.SAMPLE_GENE_LIST, hidden=True)
            custom_text_gene_list.gene_list.locked = True
            custom_text_gene_list.gene_list.save()
            # TODO: Verify gene list is good otherwise create SeqAutoRunEvent

            if custom_text_gene_list.gene_list.import_status != ImportStatus.SUCCESS:
                message = "Problem importing QC Gene List %s\n" % self.path
                message += f"Contents: {custom_gene_list_text}"
                message += f"Error: {custom_text_gene_list.gene_list.error_message}"
                logging.error(message)
                if seqauto_run:  # won't always - can be kicked off from management commands
                    SeqAutoRunEvent.objects.create(seqauto_run=seqauto_run,
                                                   file_type=SequencingFileType.QC,
                                                   message=message,
                                                   severity=LogLevel.ERROR)

            self.custom_text_gene_list = custom_text_gene_list

        try:
            # SampleFromSequencingSample is only done after VCF import, so this may not be linked yet.
            sample = self.sequencing_sample.samplefromsequencingsample.sample
            self.create_and_assign_sample_gene_list(sample)  # Also saves
        except SampleFromSequencingSample.DoesNotExist:
            self.save()

    def create_and_assign_sample_gene_list(self, sample: Sample):
        logging.info("QCGeneList: %d - create_and_assign_sample_gene_list for %s", self.pk, sample)
        self.sample_gene_list = SampleGeneList.objects.get_or_create(sample=sample,
                                                                     gene_list=self.custom_text_gene_list.gene_list)[0]
        self.save()


class QCExecSummary(SeqAutoFile, SequencingSamplePropertiesMixin):
    qc = models.ForeignKey(QC, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, null=True, on_delete=CASCADE)  # GOI - null = all genes

    percent_20x_kit = models.FloatField(null=True)
    # Coverages below are across GOIs (some of these will be null)
    percent_500x = models.FloatField(null=True)
    percent_250x = models.FloatField(null=True)
    percent_20x = models.FloatField(null=True)
    percent_10x = models.FloatField(null=True)
    mean_coverage_across_genes = models.FloatField()
    mean_coverage_across_kit = models.FloatField()
    uniformity_of_coverage = models.FloatField()
    percent_read_enrichment = models.FloatField()
    duplicated_alignable_reads = models.FloatField()
    median_insert = models.FloatField()
    ts_to_tv_ratio = models.FloatField()
    number_snps = models.IntegerField()
    snp_dbsnp_percent = models.FloatField()
    number_indels = models.IntegerField()
    indels_dbsnp_percent = models.FloatField()

    def get_coverage_columns(self):
        # TODO: Is it easier just to return non-null columns?
        HIGH_COVERAGE = ('percent_500x', 'percent_250x')
        LOW_COVERAGE = ('percent_20x', 'percent_10x')

        # If enrichment_kit is set, use that to determine columns
        enrichment_kit = self.qc.bam_file.unaligned_reads.sequencing_sample.enrichment_kit
        if enrichment_kit is not None:
            if enrichment_kit.enrichment_kit_type == EnrichmentKitType.AMPLICON:
                columns = HIGH_COVERAGE
            else:
                columns = LOW_COVERAGE

            for c in columns:  # Test
                val = getattr(self, c)
                if val is None:
                    msg = f"EnrichmentKit {enrichment_kit} columns {c} is None"
                    raise ValueError(msg)
        else:
            columns = None
            for cov in [HIGH_COVERAGE, LOW_COVERAGE]:
                all_ok = True
                for c in cov:
                    val = getattr(self, c)
                    all_ok &= val is not None
                if all_ok:
                    columns = cov
                    break
            if columns is None:
                msg = "Couldn't find the 2 non-null coverage columns in QCExecSummary: %r" % self.pk
                raise ValueError(msg)

        return columns

    @staticmethod
    def load_for_qc(_seqauto_run, qc, **kwargs):
        # We switched to loading via TSV - but had to keep the qc.path the same
        # So that it didn't muck with historical data
        exec_summary_filename = QC.get_tsv_path_from_vcf(qc.vcf_file)
        exec_summary_data = load_exec_summary(exec_summary_filename)

        # Sanity check sample names match
        # TODO: Better name or way of storing this info than "aligned_pattern"?
        params = qc.get_params()
        expected_sample_name = params['aligned_pattern']
        sample = exec_summary_data["sample"]
        if sample != expected_sample_name:
            msg = f"Exec summary file '{qc.path}' had sample name of '{sample}' while we expected '{expected_sample_name}'"
            raise ValueError(msg)

        exec_data = exec_summary_data["exec_data"]
        exec_data["qc"] = qc
        exec_summary = QCExecSummary(**exec_data)
        exec_summary.data_state = DataState.COMPLETE
        exec_summary.save()

        reference_range = exec_summary_data["reference_range"]
        if reference_range:
            reference_range["exec_summary"] = exec_summary
            esrr = ExecSummaryReferenceRange(**reference_range)
            esrr.save()

    @property
    def sequencing_sample(self):
        return self.qc.sequencing_sample

    @property
    def sample_name(self):
        return self.sequencing_sample.sample_name

    @staticmethod
    def get_sequencing_run_path():
        return "qc__bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencing_run"

    def __str__(self):
        return f"QCExecSummary ({self.data_state}) for {self.qc}"


class ExecSummaryReferenceRange(models.Model):
    exec_summary = models.OneToOneField(QCExecSummary, on_delete=CASCADE)

    percent_20x = DecimalRangeField(null=True)
    percent_10x = DecimalRangeField(null=True)
    mean_coverage_across_genes = DecimalRangeField()
    # We used to have a reference range for mean coverage, now we can have either that or a set limit (one must be set)
    mean_coverage_across_kit = DecimalRangeField(null=True)
    min_mean_coverage_across_kit = models.IntegerField(null=True)
    min_percent_20x_kit = models.IntegerField(null=True)
    uniformity_of_coverage = DecimalRangeField()
    percent_read_enrichment = DecimalRangeField()
    duplicated_alignable_reads = DecimalRangeField()
    median_insert = DecimalRangeField()
    ts_to_tv_ratio = DecimalRangeField()
    number_snps = DecimalRangeField()
    snp_dbsnp_percent = DecimalRangeField()
    number_indels = DecimalRangeField()
    indels_dbsnp_percent = DecimalRangeField()


class QCGeneCoverage(SeqAutoFile):
    qc = models.OneToOneField(QC, on_delete=CASCADE)
    gene_coverage_collection = models.OneToOneField(GeneCoverageCollection, null=True, on_delete=CASCADE)

    @staticmethod
    def get_path_from_qc(qc):
        pattern = os.path.join(qc.qc_dir_path, settings.SEQAUTO_QC_GENE_COVERAGE_PATTERN)
        return pattern % qc.get_params()

    def load_from_file(self, seqauto_run, **kwargs):
        enrichment_kit = None
        if self.qc:  # SeqAuto
            enrichment_kit = self.qc.sequencing_sample.enrichment_kit

        gene_coverage_collection = GeneCoverageCollection.objects.create(path=self.path, genome_build=self.qc.genome_build)
        self.gene_coverage_collection = gene_coverage_collection
        self.data_state = DataState.RUNNING
        self.save()

        try:
            warnings = gene_coverage_collection.load_from_file(enrichment_kit, **kwargs)
            if seqauto_run:
                for w in warnings:
                    SeqAutoRunEvent.objects.create(seqauto_run=seqauto_run,
                                                   file_type=SequencingFileType.QC,
                                                   message=w,
                                                   severity=LogLevel.WARNING)
            self.data_state = DataState.COMPLETE
        except FileNotFoundError:
            self.data_state = DataState.DELETED

        self.save()


class GoldReference(models.Model):
    """ Represents stats collected against a set of gold sequencing runs """
    enrichment_kit = models.OneToOneField(EnrichmentKit, on_delete=CASCADE)
    created = models.DateTimeField(auto_now_add=True)
    import_status = models.CharField(max_length=1, choices=ImportStatus.choices, default=ImportStatus.CREATED)
    error_exception = models.TextField(null=True)

    def __str__(self):
        import_status = self.get_import_status_display()
        name = f"{self.enrichment_kit} ({import_status})"
        if self.error_exception:
            name += " error: {self.error_exception}"
        return name


class GoldGeneCoverageCollection(models.Model):
    """ Stores what was used to calculate a gold reference (at the time)
        If you wish to delete / replace a GeneCoverageCollection here, you must delete the
        old gold first (PROTECT) to stop Stored gold reference and current gold runs
        getting out of sync """
    SEQUENCING_SAMPLE_PATH = "gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample"

    gold_reference = models.ForeignKey(GoldReference, on_delete=CASCADE)
    gene_coverage_collection = models.ForeignKey(GeneCoverageCollection, on_delete=PROTECT)

    @property
    def sequencing_sample(self):
        try:
            ss = self.gene_coverage_collection.qcgenecoverage.qc.bam_file.unaligned_reads.sequencing_sample
        except:
            ss = None
        return ss

    @property
    def sequencing_run(self):
        ss = self.sequencing_sample
        if ss:
            return self.sequencing_sample.sample_sheet.sequencing_run
        return None


class GoldCoverageSummary(models.Model):
    gold_reference = models.ForeignKey(GoldReference, on_delete=CASCADE)
    gene_symbol = models.ForeignKey(GeneSymbol, null=True, on_delete=CASCADE)
    transcript = models.ForeignKey(Transcript, null=True, on_delete=CASCADE)
    original_gene_symbol = models.TextField()
    original_transcript_id = models.TextField()
    mean = models.FloatField()
    standard_error = models.FloatField()
    min_mean = models.FloatField()
    depth_20x_5th_percentile = models.FloatField()
    depth_10x_5th_percentile = models.FloatField()
    depth_mean_5th_percentile = models.FloatField()
    depth_mean_95th_percentile = models.FloatField()

    @staticmethod
    def filter_for_kit_and_gene_symbols(enrichment_kit, gene_symbols):
        return GoldCoverageSummary.objects.filter(gene_symbol__in=gene_symbols,
                                                  gold_reference__enrichment_kit=enrichment_kit)


class QCType(models.Model):
    name = models.TextField(primary_key=True)
    qc_object_path = models.TextField()  # From QC object
    total_field = models.TextField(null=True)  # If present, be able to show % vs this field

    def __str__(self):
        return self.name


class QCColumn(models.Model):
    qc_type = models.ForeignKey(QCType, on_delete=CASCADE)
    name = models.TextField()
    field = models.TextField()  # field from 'qc_type.qc_object_path'

    def __str__(self):
        return self.name


class JobScript(models.Model):
    # Attach the samplesheet/bam/vcf/qc to this script so delete will cascade through
    FIELDS = {
        SequencingFileType.SAMPLE_SHEET: 'sample_sheet',  # TODO: Rename to basecalling?
        SequencingFileType.ILLUMINA_FLOWCELL_QC: 'illumina_flowcell_qc',
        SequencingFileType.FASTQC: 'fastqc',
        SequencingFileType.BAM: 'bam_file',
        SequencingFileType.FLAGSTATS: 'flagstats',
        SequencingFileType.VCF: 'vcf_file',
        SequencingFileType.COMBINED_VCF: 'combined_vcf_file',
        SequencingFileType.QC: 'qc',
        SequencingFileType.DATA_MIGRATION: 'sample_sheet',
    }

    seqauto_run = models.ForeignKey(SeqAutoRun, on_delete=CASCADE)
    path = models.TextField()  # script
    out_file = models.TextField(null=True)
    file_type = models.CharField(max_length=1, choices=SequencingFileType.choices)
    sample_sheet = models.ForeignKey(SampleSheet, null=True, on_delete=CASCADE)
    illumina_flowcell_qc = models.ForeignKey(IlluminaFlowcellQC, null=True, on_delete=CASCADE)
    fastqc = models.ForeignKey(FastQC, null=True, on_delete=CASCADE)
    bam_file = models.ForeignKey(BamFile, null=True, on_delete=CASCADE)
    flagstats = models.ForeignKey(Flagstats, null=True, on_delete=CASCADE)
    vcf_file = models.ForeignKey(VCFFile, null=True, on_delete=CASCADE)
    combined_vcf_file = models.ForeignKey(SampleSheetCombinedVCFFile, null=True, on_delete=CASCADE)
    qc = models.ForeignKey(QC, null=True, on_delete=CASCADE)
    job_id = models.TextField(null=True)
    job_status = models.CharField(max_length=1, choices=JobScriptStatus.choices, default=JobScriptStatus.CREATED)
    return_code = models.IntegerField(null=True)

    def get_record(self):
        f = dict(JobScript.FIELDS)[self.file_type]
        return getattr(self, f)

    def get_variable_name(self):
        return f"{self.file_type}_{self.pk}"

    def job_complete(self):
        self.job_status = JobScriptStatus.FINISHED
        self.save()
        error_exception = None
        record = self.get_record()

        # TODO: Hmmm, maybe this is a bit dodgy. Perhaps make a new DataState which is
        # Unprocessed, and during save() check if that, then do loading and set to complete?
        if self.return_code:
            logging.error("JobScript %s returned %d", self, self.return_code)
            data_state = DataState.ERROR
        else:
            if hasattr(record, 'load_from_file'):
                try:
                    record.load_from_file(self.seqauto_run)
                    data_state = DataState.COMPLETE
                except:
                    error_exception = get_traceback()
                    data_state = DataState.ERROR
            else:
                data_state = DataState.COMPLETE

        record.error_exception = error_exception
        record.data_state = data_state
        record.save()

    def __str__(self):
        record = self.get_record()
        record_pk = record.pk if record else 'N/A'
        return "%s: %r" % (self.file_type, record_pk)


class SequencingInfo(models.Model):
    lab_project = models.ForeignKey(LabProject, on_delete=CASCADE)
    doi = models.TextField(blank=True, null=True)
    paper_name = models.TextField(blank=True, null=True)
    year_published = models.IntegerField(null=True)
    enrichment_kit = models.ForeignKey(EnrichmentKit, on_delete=CASCADE, blank=True, null=True)
    sequencer = models.ForeignKey(Sequencer, on_delete=CASCADE, blank=True, null=True)
    seq_details = models.TextField(blank=True, null=True)
    file_type = models.TextField(blank=True, null=True)
    file_count = models.IntegerField(blank=True, default=0)


@receiver(post_delete, sender=JobScript)
def post_delete_job_script(sender, instance, *args, **kwargs):
    try:
        os.remove(instance.path)
    except:
        pass


def get_samples_by_sequencing_sample(sample_sheet, vcf):
    sequencing_samples_by_name = sample_sheet.get_sequencing_samples_by_name()

    def clean_sample_name(s):
        return s.upper().replace("-", "_")

    samples_by_sequencing_sample = {}
    for sample in vcf.sample_set.all():
        # Do a startswith match rather than hash lookup as it's less strict (diff naming conventions etc)
        for sequencing_sample_name, sequencing_sample in sequencing_samples_by_name.items():
            cleaned_sample_name = clean_sample_name(sample.name)
            cleaned_sequencing_sample_name = clean_sample_name(sequencing_sample_name)
            if cleaned_sample_name.startswith(cleaned_sequencing_sample_name):
                samples_by_sequencing_sample[sequencing_sample] = sample
    return samples_by_sequencing_sample


def get_20x_gene_coverage(gene_symbol, min_coverage=100):
    unique_gene_coverage_qs = gene_symbol.genecoveragecanonicaltranscript_set.filter(gene_coverage_collection__qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencingruncurrentsamplesheet__isnull=False)
    unique_gene_coverage_qs = unique_gene_coverage_qs.distinct("gene_coverage_collection")
    return unique_gene_coverage_qs.filter(percent_20x__gte=min_coverage).count()
