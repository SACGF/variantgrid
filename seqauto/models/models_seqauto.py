import logging
import os
import pathlib
import re
import shutil
from datetime import datetime
from functools import cached_property
from typing import Optional

from django.conf import settings
from django.contrib import messages
from django.contrib.auth.models import User, Group
from django.contrib.postgres.fields import DecimalRangeField
from django.core.cache import cache
from django.db import models
from django.db.models import Value, When, Case, IntegerField, Max
from django.db.models.deletion import SET_NULL, CASCADE, PROTECT
from django.db.models.signals import post_delete, pre_delete
from django.dispatch.dispatcher import receiver
from django.urls.base import reverse
from django.utils.timezone import make_aware
from django_extensions.db.models import TimeStampedModel

from genes.models import GeneListCategory, CustomTextGeneList, GeneList, GeneCoverageCollection, \
    Transcript, GeneSymbol, SampleGeneList, TranscriptVersion, GeneCoverageCanonicalTranscript, ActiveSampleGeneList
from library.constants import DAY_SECS
from library.enums.log_level import LogLevel
from library.genomics.vcf_utils import get_variant_caller_and_version_from_vcf
from library.log_utils import get_traceback, log_traceback
from library.preview_request import PreviewModelMixin
from library.utils import sorted_nicely
from library.utils.file_utils import name_from_filename, remove_gz_if_exists
from patients.models import FakeData, Patient
from seqauto.illumina import illuminate_report
from seqauto.illumina.illumina_sequencers import SEQUENCING_RUN_REGEX
from seqauto.models.models_enums import DataGeneration, SequencerRead, PairedEnd, \
    SequencingFileType, JobScriptStatus, SeqAutoRunStatus
from seqauto.models.models_sequencing import Sequencer, EnrichmentKit, Experiment
from seqauto.models.models_software import Aligner, VariantCaller
from seqauto.qc.exec_summary import load_exec_summary
from seqauto.qc.fastqc_parser import read_fastqc_data
from seqauto.qc.flag_stats import load_flagstats
from seqauto.qc.qc_utils import meta_data_file
from seqauto.signals.signals_list import sequencing_run_sample_sheet_created_signal
from snpdb.models import VCF, Sample, GenomeBuild, DataState, InheritanceManager, Wiki
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

    class Meta:
        permissions = (('seqauto_scan_initiate', 'SeqAuto scan initiate'),)

    def save(self, *args, **kwargs):
        self.status = self.get_status()
        super().save(*args, **kwargs)

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


class SeqAutoRecord(TimeStampedModel):
    """ Base class for everything below """
    objects = InheritanceManager()
    path = models.TextField()
    # Everything is derived from a SequencingRun so just keep this around to simplify traversing models
    # Need to keep this nullable so we can make a SequencingRun before assigning FK to itself
    sequencing_run = models.ForeignKey("SequencingRun", null=True, on_delete=CASCADE)
    # Stores stat st_mtime - time of last modification - only used for classes that can reload
    file_last_modified = models.FloatField(default=0.0)
    hash = models.TextField(blank=True)  # Not used for everything
    is_valid = models.BooleanField(default=False)  # Set in save
    # data_state was used to create 'expected' objects ie vcfs for bam files
    # that was then set based on whether the file turned up. If it disappeared it would be set to DELETED
    # But with API - we assume anything sent to us is COMPLETED
    # We will probably remove this field in the future as we go API only
    data_state = models.CharField(max_length=1, choices=DataState.choices)

    def save(self, *args, **kwargs):
        self.validate()
        self.is_valid = not SeqAutoMessage.objects.filter(record=self, open=True, severity=LogLevel.ERROR).exists()
        super().save(*args, **kwargs)

    def validate(self) -> bool:
        """ Creates SeqAutoMessage of severity=ERROR - or closes any that no longer apply """
        return True

    @property
    def last_modified_datetime(self):
        return make_aware(datetime.fromtimestamp(self.file_last_modified))

    def _close_messages_with_code(self, *args):
        SeqAutoMessage.objects.filter(record=self, code__in=args).update(open=False)

    @staticmethod
    def get_file_last_modified(filename):
        path = pathlib.Path(filename)
        return path.stat().st_mtime

    def add_messages(self, request):
        qs = SeqAutoMessage.objects.filter(record=self, open=True).annotate(
            priority=Case(
                When(severity=LogLevel.ERROR, then=Value(1)),
                When(severity=LogLevel.WARNING, then=Value(2)),
                When(severity=LogLevel.INFO, then=Value(3)),
                When(severity=LogLevel.DEBUG, then=Value(4)),
                output_field=IntegerField(),
            )
        ).order_by("priority")
        for sa_message in qs:
            messages.add_message(request, sa_message.level, sa_message.message)


class SeqAutoMessage(TimeStampedModel):
    # TODO: Remember we can bump to latest SeqAutoRun as we go
    # At least one of seqauto_run/record should be set
    seqauto_run = models.ForeignKey(SeqAutoRun, null=True, on_delete=CASCADE)
    record = models.ForeignKey(SeqAutoRecord, null=True, on_delete=CASCADE)
    severity = models.CharField(max_length=1, choices=LogLevel.CHOICES)
    code = models.TextField(null=True)  # A code you can use to close messages after resolving
    message = models.TextField()
    open = models.BooleanField(default=True)

    @property
    def level(self):
        """ Convert from LogLevel to Django Messages constants """
        return messages.constants.DEFAULT_LEVELS[self.get_severity_display()]

    def __str__(self):
        record = SeqAutoRecord.objects.get_subclass(pk=self.record)
        return f"{self.seqauto_run} {self.get_severity_display()} {record}: {self.message}"


class SequencingRun(PreviewModelMixin, SeqAutoRecord):
    """ Represents a flowcell (or other technology with multiple sequencing samples) """
    name = models.TextField(primary_key=True)
    date = models.DateField(null=True)  # e.g. from the sequencing name - used to sort
    sequencer = models.ForeignKey(Sequencer, on_delete=CASCADE)
    gold_standard = models.BooleanField(default=False)
    bad = models.BooleanField(default=False)
    hidden = models.BooleanField(default=False)
    legacy = models.BooleanField(default=False)  # Don't update it with scans (eg to say file missing)
    experiment = models.ForeignKey(Experiment, null=True, blank=True, on_delete=SET_NULL)
    # Sequencing Run can be all one enrichment_kit, or SequencingSample can have own enrichment_kits
    enrichment_kit = models.ForeignKey(EnrichmentKit, null=True, blank=True, on_delete=CASCADE)
    has_basecalls = models.BooleanField(default=False)
    has_interop = models.BooleanField(default=False)  # Quality, Index and Tile
    fake_data = models.ForeignKey(FakeData, null=True, blank=True, on_delete=CASCADE)

    @classmethod
    def preview_icon(cls) -> str:
        return "fa-solid fa-person-running"

    @classmethod
    def preview_enabled(cls) -> bool:
        return settings.SEQAUTO_ENABLED

    def _validate(self):
        sample_sheet_changed_code = "sample_sheet_changed"

        if self.is_data_out_of_date_from_current_sample_sheet:
            SeqAutoMessage.objects.update_or_create(record=self, severity=LogLevel.ERROR,
                                                    code=sample_sheet_changed_code,
                                                    message="SampleSheet has changed, please confirm in 'Admin' tab",
                                                    defaults={"open": True})
        else:
            self._close_messages_with_code(sample_sheet_changed_code)

        sample_sheet_data_state = "sample_sheet_data_state"
        sample_sheet_missing = "sample_sheet_missing"
        sample_sheet_qc_exception = "sample_sheet_qc_exception"
        try:
            illumina_qc = self.sequencingruncurrentsamplesheet.sample_sheet.illuminaflowcellqc
            self._close_messages_with_code(sample_sheet_missing, sample_sheet_qc_exception)

            if illumina_qc.data_state != DataState.COMPLETE:
                SeqAutoMessage.objects.update_or_create(record=self, severity=LogLevel.ERROR,
                                                        code=sample_sheet_data_state,
                                                        message=f"QC data state={illumina_qc.get_data_state_display()}",
                                                        defaults={"open": True})
            else:
                self._close_messages_with_code(sample_sheet_data_state)

        except SequencingRunCurrentSampleSheet.DoesNotExist:
            SeqAutoMessage.objects.update_or_create(record=self, severity=LogLevel.ERROR,
                                                    code=sample_sheet_missing,
                                                    message="No Current SampleSheet set",
                                                    defaults={"open": True})
        except Exception as e:
            SeqAutoMessage.objects.update_or_create(record=self, severity=LogLevel.ERROR,
                                                    sample_sheet_qc_exception=sample_sheet_qc_exception,
                                                    message=f"QC not loaded: {e}",
                                                    defaults={"open": True})

    def get_current_sample_sheet(self):
        try:
            return self.sequencingruncurrentsamplesheet.sample_sheet
        except:
            logging.info("Can't find current sample sheet for %s", self.path)
            raise

    @staticmethod
    def get_name_validation_errors(name: str) -> Optional[str]:
        """ Returns a string of issues with the name """
        errors = None
        if settings.SEQAUTO_SEQUENCING_RUN_VALIDATE_ILLUMINA_FORMULA:
            if not re.search(SEQUENCING_RUN_REGEX, name):
                errors = f"Name does not match Illumina regex: '{SEQUENCING_RUN_REGEX}'"
        return errors

    @staticmethod
    def get_original_illumina_sequencing_run(modified_sequencing_run):
        # TAU rename the sequencing run dir with enrichment kit at the end - need to clean it
        original_sequencing_run = modified_sequencing_run
        if m := re.search(SEQUENCING_RUN_REGEX, modified_sequencing_run):
            original_sequencing_run = m.group(0)
        return original_sequencing_run

    @staticmethod
    def get_date_from_name(name) -> Optional[datetime.date]:
        date = None
        if m := re.match(r".*_?([12]\d{5})_", name):
            date_str = m.group(1)
            dt = datetime.strptime(date_str, "%y%m%d")
            date = make_aware(dt).date()
        return date

    def get_params(self):
        """ This allows chaining down names etc - in case a level changes it, will cascade down """
        params = {
            "sequencing_run": self.name,
            "sequencing_run_dir": self.path,
            "original_sequencing_run": self.get_original_illumina_sequencing_run(self.name),
        }
        if self.enrichment_kit:
            params["enrichment_kit"] = self.enrichment_kit.name

        if self.experiment:
            params["experiment"] = self.experiment.name

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
        interop_dir = os.path.join(self.path, settings.SEQAUTO_SEQUENCING_RUN_INTEROP_SUB_DIR)
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

    @property
    def is_data_out_of_date_from_current_sample_sheet(self):
        try:
            current_sample_sheet = self.get_current_sample_sheet()
        except SequencingRunCurrentSampleSheet.DoesNotExist:
            return False

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


class SampleSheet(SeqAutoRecord):
    SAMPLE_SHEET = "SampleSheet.csv"

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
        date = "TODO"
        return f"{self.hash}/{date}"

    def set_as_current_sample_sheet(self, sequencing_run, created=False, seqauto_run=None):
        if created:
            sequencing_run_sample_sheet_created_signal.send(sender=os.path.basename(__file__),
                                                            sample_sheet=self)

        # Make sure SequencingRunCurrentSampleSheet is set to what we found on disk
        try:
            # Update existing
            current_ss = sequencing_run.sequencingruncurrentsamplesheet
            on_disk_not_current = current_ss.sample_sheet != self
            if on_disk_not_current:
                from seqauto.sequencing_files.create_resource_models import current_sample_sheet_changed
                current_sample_sheet_changed(seqauto_run, current_ss, self)

        except SequencingRunCurrentSampleSheet.DoesNotExist:
            # Create new
            current_ss = SequencingRunCurrentSampleSheet.objects.create(sequencing_run=sequencing_run,
                                                                        sample_sheet=self)
            logging.info("Created new SequencingRunCurrentSampleSheet: %s", current_ss)

    def __str__(self):
        return self.path


class SequencingRunCurrentSampleSheet(models.Model):
    """ This is a way to "walk" via relations from sequencing run to the latest sample sheet """
    sequencing_run = models.OneToOneField(SequencingRun, on_delete=CASCADE)
    sample_sheet = models.OneToOneField(SampleSheet, on_delete=CASCADE)

    def __str__(self):
        return f"{self.sequencing_run}, {self.sample_sheet}"


class SequencingSample(models.Model):
    """ Represents a row in a SampleSheet.csv

        As it's not a file, not a SeqAutoRecord
     """
    sample_sheet = models.ForeignKey(SampleSheet, on_delete=CASCADE)
    sample_id = models.TextField()
    # sample_name is used to name files. In MiSeq/NextSeq samplesheet you can add names.
    # For Hiseq and if left empty on MiSeq this will be sample_id
    sample_name = models.TextField(null=True)
    sample_project = models.TextField(null=True)
    sample_number = models.IntegerField()  # Row from sample sheet
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

    @cached_property
    def patient(self) -> Optional[Patient]:
        patients_set = set()
        for sfss in self.samplefromsequencingsample_set.all():
            if sfss.sample.patient:
                patients_set.add(sfss.sample.patient)

        patient = None
        if patients_set:
            if len(patients_set) > 1:
                raise ValueError(f"{self} is linked multiple samples with different patients!")
            patient = patients_set.pop()
        return patient

    def __str__(self):
        return self.sample_id


class SequencingSampleData(models.Model):
    """ key/values set from settings.SEQAUTO_SAMPLE_SHEET_EXTRA_COLUMNS """
    sequencing_sample = models.ForeignKey(SequencingSample, on_delete=CASCADE)
    column = models.TextField()
    value = models.TextField(null=True)


class SampleFromSequencingSample(models.Model):
    sample = models.OneToOneField(Sample, on_delete=CASCADE)
    sequencing_sample = models.ForeignKey(SequencingSample, on_delete=CASCADE)

    @property
    def sequencing_run(self):
        return self.sequencing_sample.sample_sheet.sequencing_run


class VCFFromSequencingRun(models.Model):
    """ This object exists so it's easy to build VCF links on the sequencing list grid """
    vcf = models.OneToOneField(VCF, on_delete=CASCADE)
    sequencing_run = models.ForeignKey(SequencingRun, on_delete=CASCADE)
    variant_caller = models.ForeignKey(VariantCaller, null=True, on_delete=SET_NULL)

    class Meta:
        unique_together = ("sequencing_run", "variant_caller")

    def get_variant_caller(self):
        """ Because we can't do a schema change for VG3, we need to basically always store variant_caller as None
            for single sample VCFs """
        if self.variant_caller:
            return self.variant_caller
        try:
            return self.vcf.uploadedvcf.backendvcf.variant_caller
        except:
            pass
        return None


class IlluminaFlowcellQC(SeqAutoRecord):
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
        illuminate_dir = settings.SEQAUTO_ILLUMINATE_QC_DIR_PATTERN % sequencing_run.get_params()  # this may fail if
        # sa path signals haven't handled enrichment kit
        return os.path.join(illuminate_dir, "illuminate_report.txt")

    @staticmethod
    def get_sequencing_run_path():
        return 'sample_sheet__sequencing_run'

    def __str__(self):
        return self.path


class ReadQ30(models.Model):
    illumina_flowcell_qc = models.ForeignKey(IlluminaFlowcellQC, on_delete=CASCADE)
    sequencer_read_id = models.IntegerField()  # e.g. HiSeq= [R1,Index,R2], NextSeq/MiSeq=[R1,Index1,Index2,R2]
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


class Fastq(SeqAutoRecord):
    sequencing_sample = models.ForeignKey(SequencingSample, on_delete=CASCADE)
    name = models.TextField()  # from path
    read = models.CharField(max_length=2, choices=PairedEnd.choices)

    def get_common_pair_name_and_read(self):
        """ Need a way to group R1/R2 FastQs together (using globally unique key - ie shared path prefix)  """
        regex_pattern = "(.+)_(R[12])(_001)?.fastq.gz"
        m = re.match(regex_pattern, self.path)
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
        """ The 1st value returned must be what will be generated given current pipelines
            (ie old naming conventions come after) """

        # Old code had FastQ names auto-generated by Illumina basecalling:
        # sequencer_model = sequencing_sample.sample_sheet.sequencing_run.sequencer.sequencer_model
        # patterns = ["%(full_sample_name)s_L%(lane)03d_R%(read)d_001.fastq.gz"]
        #
        # if sequencer_model.data_naming_convention == DataGeneration.HISEQ:
        #     patterns.append("%(sample_id)s_R%(read)d.fastq.gz")  # SACGF style naming convention
        # elif sequencer_model.data_naming_convention == DataGeneration.MISEQ:
        #     patterns.extend(["%(sample_name)s_S%(sample_number)d_R%(read)d_001.fastq.gz",
        #                      "%(sample_name_underscores)s_S%(sample_number)d_L%(lane)03d_R%(read)d_001.fastq.gz",
        #                      "%(sample_name_underscores)s_S%(sample_number)d_R%(read)d_001.fastq.gz"])
        # else:
        #     msg = f"Unknown sequencer_model.data_naming_convention '{sequencer_model.data_naming_convention}'"
        #     raise ValueError(msg)

        # New Diagnostic pipeline all FastQs have this simple format
        patterns = ["%(sample_id)s_R%(read)d.fastq.gz"]

        params = sequencing_sample.get_params()
        r1_params = {"read": 1}
        r1_params.update(params)
        r2_params = {"read": 2}
        r2_params.update(params)

        unaligned_dir_patterns = []
        if sequencing_sample.sample_project:
            fastq_subdirs = ["%(sequencing_run_dir)s/Data/Intensities/BaseCalls/%(sample_project)s",
                             "%(sequencing_run_dir)s/Data/Intensities/BaseCalls/%(sample_project)s/%(sample_id)s"]
            unaligned_dir_patterns.extend(fastq_subdirs)
        unaligned_dir_patterns.append("%(sequencing_run_dir)s/Data/Intensities/BaseCalls")
        unaligned_dir_patterns.append(settings.SEQAUTO_FASTQ_DIR_PATTERN)

        pair_paths = []
        for unaligned_dir_pattern in unaligned_dir_patterns:
            try:
                pattern = unaligned_dir_pattern
                fastq_dir = pattern % params

                for pattern in patterns:
                    name1 = pattern % r1_params
                    name2 = pattern % r2_params

                    fastq_r1 = os.path.abspath(os.path.join(fastq_dir, name1))
                    fastq_r2 = os.path.abspath(os.path.join(fastq_dir, name2))
                    pair_paths.append((fastq_r1, fastq_r2))
            except KeyError as ke:
                logging.error("Cannot populate old style string formatting: '%s', params: %s: %s",
                              pattern, params, ke)
                raise

        return pair_paths

    def __str__(self):
        return f"FastQ {self.name} ({self.read}) from sample {self.sequencing_sample}"


class FastQC(SeqAutoRecord):
    fastq = models.OneToOneField(Fastq, on_delete=CASCADE)
    total_sequences = models.IntegerField(null=True)
    filtered_sequences = models.IntegerField(null=True)
    gc = models.IntegerField(null=True)

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
        fastqc_dir = f"{name_from_filename(name)}_fastqc"
        path = os.path.join(base_dir, "FastQC", fastqc_dir, "fastqc_data.txt")
        return os.path.abspath(path)

    def __str__(self):
        return f"FastQC ({self.get_data_state_display()}) for {self.fastq}"


class UnalignedReads(models.Model):
    """ Not a file just a way of keeping paired end fastqs together """
    sequencing_sample = models.ForeignKey(SequencingSample, on_delete=CASCADE)
    fastq_r1 = models.ForeignKey(Fastq, related_name='fastq_r1', on_delete=CASCADE)
    fastq_r2 = models.ForeignKey(Fastq, null=True, related_name='fastq_r2', on_delete=CASCADE)

    @property
    def sequencing_run(self):
        return self.sequencing_sample.sample_sheet.sequencing_run

    def get_params(self):
        fastq_params = self.fastq_r1.sequencing_sample.get_params()
        sequencing_run = self.sequencing_sample.sample_sheet.sequencing_run
        data_naming_convention = sequencing_run.sequencer.sequencer_model.data_naming_convention
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


class BamFile(SeqAutoRecord):
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
        pattern = os.path.join(settings.SEQAUTO_BAM_DIR_PATTERN, settings.SEQAUTO_BAM_PATTERN)
        try:
            filename = pattern % params
        except KeyError as ke:
            if "enrichment_kit" in ke.args:
                logging.error("'enrichment_kit' not set, this is usually done via signal handlers in a custom app")
            missing = ', '.join(ke.args)
            logging.error("%s missing: %s. params=%s", unaligned_reads, missing, params)
            raise
        return os.path.abspath(filename)

    @staticmethod
    def get_aligner_from_bam_file(bam_path):
        # TODO: Do properly
        aligner, _ = Aligner.objects.get_or_create(name='fake_aligner', version='0.666')
        return aligner

    def __str__(self):
        return f"BAM: {self.unaligned_reads.sequencing_sample} ({self.get_data_state_display()})"


class Flagstats(SeqAutoRecord):
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


class VCFFile(SeqAutoRecord):
    """ VCFs from the file system """
    bam_file = models.ForeignKey(BamFile, on_delete=CASCADE)
    variant_caller = models.ForeignKey(VariantCaller, on_delete=CASCADE)

    def get_params(self):
        return self.bam_file.get_params()

    @property
    def name(self):
        return name_from_filename(self.path)

    @property
    def upload_pipeline(self):
        try:
            up = self.backendvcf.uploaded_vcf.upload_pipeline
        except:
            up = None
        return up

    @property
    def sample_sheet(self) -> SampleSheet:
        return self.bam_file.unaligned_reads.sequencing_sample.sample_sheet

    def load_from_file(self, seqauto_run, **kwargs):
        if not settings.SEQAUTO_IMPORT_VCF:
            raise DontAutoLoadException()

        import_filesystem_vcf(self.path, self.name, ImportSource.SEQAUTO)

    @staticmethod
    def get_path_from_bam(bam):
        vcf_pattern = None
        if enrichment_kit := bam.sequencing_sample.enrichment_kit:
            vcf_pattern = settings.SEQAUTO_VCF_PATTERNS_FOR_KIT.get(enrichment_kit.name)
        if vcf_pattern is None:
            vcf_pattern = settings.SEQAUTO_VCF_PATTERNS_FOR_KIT["default"]
        pattern = os.path.join(settings.SEQAUTO_VCF_DIR_PATTERN, vcf_pattern)
        vcf_path = pattern % bam.get_params()
        return os.path.abspath(vcf_path)

    def __str__(self):
        return f"VCF {self.name} ({self.get_data_state_display()})"


class SampleSheetCombinedVCFFile(SeqAutoRecord):
    sample_sheet = models.ForeignKey(SampleSheet, on_delete=CASCADE)
    variant_caller = models.ForeignKey(VariantCaller, on_delete=CASCADE)

    def get_params(self):
        return self.sample_sheet.get_params()

    @property
    def name(self):
        """ What the VCF gets named """
        return os.path.basename(self.path)

    @property
    def upload_pipeline(self):
        try:
            up = self.backendvcf.uploaded_vcf.upload_pipeline
        except:
            up = None
        return up

    @cached_property
    def vcf(self) -> Optional[VCF]:
        try:
            VCF.objects.get(uploadedvcf__uploaded_file__path=self.path)
        except VCF.DoesNotExist:
            return None

    def load_from_file(self, seqauto_run, **kwargs):
        if not settings.SEQAUTO_IMPORT_COMBO_VCF:
            raise DontAutoLoadException()

        if self.vcf:
            # Already exists! Make sure it's linked but don't reload
            try:
                from upload.models import UploadedVCF
                from upload.vcf.vcf_import import create_backend_vcf_links

                uploaded_vcf = UploadedVCF.objects.get(uploaded_file__path=self.path)
                create_backend_vcf_links(uploaded_vcf)
            except:
                log_traceback()

            # Re-linking SequencingRun / backend VCF will be done manually in view_sequencing_run
            # "Assign data to most recent sample sheet" button
            logging.info("Skipping loading of %s as vcf named already exists!", self.path)
            raise DontAutoLoadException()

        import_filesystem_vcf(self.path, self.name, ImportSource.SEQAUTO)

    def needs_to_be_linked(self):
        try:
            _ = self.backendvcf  # if no exception we can link it
            try:
                _ = self.backendvcf.vcf.vcffromsequencingrun
            except VCFFromSequencingRun.DoesNotExist:
                return True
        except:
            pass
        return False

    @staticmethod
    def get_paths_from_sample_sheet(sample_sheet) -> list[str]:
        enrichment_kit = sample_sheet.sequencing_run.enrichment_kit
        pattern_list = None
        if enrichment_kit:
            kit_name = enrichment_kit.name
            if "ffpe" not in kit_name.lower() and "_FFPE_" in sample_sheet.sequencing_run.name:
                kit_name += "_ffpe"
            pattern_list = settings.SEQAUTO_COMBINED_VCF_PATTERNS_FOR_KIT.get(kit_name)

        if pattern_list is None:
            pattern_list = settings.SEQAUTO_COMBINED_VCF_PATTERNS_FOR_KIT["default"]

        paths = []
        for pattern in pattern_list:
            combo_path = os.path.join(settings.SEQAUTO_VCF_DIR_PATTERN, pattern) % sample_sheet.get_params()
            paths.append(os.path.abspath(combo_path))
        return paths

    def __str__(self):
        num_samples = self.sample_sheet.sequencingsample_set.count()
        return f"{self.variant_caller} ComboVCF ({self.pk}) for {self.sequencing_run} ({num_samples} samples)"


class QC(SeqAutoRecord):
    """ This holds together all of the other QC objects """
    bam_file = models.ForeignKey(BamFile, on_delete=CASCADE)
    vcf_file = models.ForeignKey(VCFFile, on_delete=CASCADE)

    @property
    def genome_build(self):
        gb = GenomeBuild.objects.filter(vcf__sample__samplefromsequencingsample__sequencing_sample=self.sequencing_sample).first()
        if gb is None:
            logging.warning("%s: requested genome build, but don't know (no VCFs linked)", self)
            gb = GenomeBuild.legacy_build()
        return gb

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


class QCGeneList(SeqAutoRecord):
    """ This represents a text file containing genes which will be used for initial pass and QC filters

        The VCF may not be loaded yet, so we'll just store gene list and link it later

        The reason we have both a sample_gene_list and a custom_text_gene_list is because we wanted to
        represent the text from a file on disk.

        I think we probably could have gotten around that as
        SeqAutoRecord contains a hash - could maybe remove that and just use gene list
    """
    qc = models.ForeignKey(QC, on_delete=CASCADE)
    custom_text_gene_list = models.OneToOneField(CustomTextGeneList, null=True, on_delete=SET_NULL)
    sample_gene_list = models.ForeignKey(SampleGeneList, null=True, on_delete=SET_NULL)

    @property
    def sequencing_sample(self):
        return self.qc.sequencing_sample

    @staticmethod
    def get_path_from_qc(qc):
        pattern = os.path.join(settings.SEQAUTO_GOI_DIR_PATTERN, settings.SEQAUTO_GOI_LIST_PATTERN)
        return pattern % qc.get_params()

    @staticmethod
    def create_gene_list(custom_gene_list_text, sequencing_sample):
        from genes.custom_text_gene_list import create_custom_text_gene_list

        custom_text_gene_list = CustomTextGeneList(name=f"QC GeneList for {sequencing_sample.sample_name}",
                                                   text=custom_gene_list_text)
        custom_text_gene_list.save()
        seqauto_user = get_seqauto_user()
        create_custom_text_gene_list(custom_text_gene_list, seqauto_user, GeneListCategory.SAMPLE_GENE_LIST,
                                     hidden=True)
        custom_text_gene_list.gene_list.locked = True
        custom_text_gene_list.gene_list.save()
        return custom_text_gene_list

    def load_from_file(self, seqauto_run, **kwargs):
        with open(self.path, encoding="utf-8") as f:
            custom_gene_list_text = f.read()
            custom_text_gene_list = self.create_gene_list(custom_gene_list_text, self.sequencing_sample)
            if custom_text_gene_list.gene_list.import_status != ImportStatus.SUCCESS:
                message = f"Problem importing QC Gene List {self.path}\n"
                message += f"Contents: {custom_gene_list_text}"
                message += f"Error: {custom_text_gene_list.gene_list.error_message}"
                SeqAutoMessage.objects.create(seqauto_run=seqauto_run,
                                              record=self,
                                              message=message,
                                              severity=LogLevel.ERROR)

            self.custom_text_gene_list = custom_text_gene_list
            self.save()

        self.link_samples_if_exist()

    def link_samples_if_exist(self, force_active=False):
        # SampleFromSequencingSample is only done after VCF import, so this may not be linked yet.
        for sfss in self.sequencing_sample.samplefromsequencingsample_set.all():
            sample = sfss.sample
            self.create_and_assign_sample_gene_list(sample, force_active=force_active)  # Also saves

    def create_and_assign_sample_gene_list(self, sample: Sample, force_active=False):
        logging.info("QCGeneList: %d - create_and_assign_sample_gene_list for %s", self.pk, sample)
        # On create, 'sample_gene_list_created' sets to active gene list
        self.sample_gene_list = SampleGeneList.objects.get_or_create(sample=sample,
                                                                     gene_list=self.custom_text_gene_list.gene_list)[0]
        self.save()

        if force_active:
            ActiveSampleGeneList.objects.update_or_create(sample=sample,
                                                          defaults={'sample_gene_list': self.sample_gene_list})


class QCExecSummary(SeqAutoRecord):
    qc = models.ForeignKey(QC, on_delete=CASCADE)
    gene_list = models.ForeignKey(GeneList, null=True, on_delete=CASCADE)  # GOI - null = all genes

    deduplicated_reads = models.IntegerField(null=True)
    indels_dbsnp_percent = models.FloatField(null=True)
    mean_coverage_across_genes = models.FloatField()
    mean_coverage_across_kit = models.FloatField()
    median_insert = models.FloatField()
    number_indels = models.IntegerField(null=True)
    number_snps = models.IntegerField(null=True)
    percent_10x_goi = models.FloatField(null=True)
    percent_20x_goi = models.FloatField(null=True)
    percent_20x_kit = models.FloatField(null=True)
    percent_100x_goi = models.FloatField(null=True)
    percent_100x_kit = models.FloatField(null=True)
    percent_250x_goi = models.FloatField(null=True)
    percent_250x_kit = models.FloatField(null=True)
    percent_500x_goi = models.FloatField(null=True)
    percent_500x_kit = models.FloatField(null=True)
    percent_error_rate = models.FloatField(null=True)
    percent_map_to_diff_chr = models.FloatField(null=True)
    percent_read_enrichment = models.FloatField()
    percent_reads = models.FloatField(null=True)
    percent_softclip = models.FloatField(null=True)
    percent_duplication = models.FloatField()
    reads = models.IntegerField(null=True)
    sample_id_lod = models.FloatField(null=True)
    sex_match = models.TextField(null=True)
    snp_dbsnp_percent = models.FloatField(null=True)
    ts_to_tv_ratio = models.FloatField(null=True)
    uniformity_of_coverage = models.FloatField()

    def get_coverage_columns(self) -> list[str]:
        COVERAGE_COLUMNS = [
            "percent_10x_goi",
            "percent_20x_goi",
            "percent_20x_kit",
            "percent_100x_goi",
            "percent_100x_kit",
            "percent_250x_goi",
            "percent_250x_kit",
            "percent_500x_goi",
            "percent_500x_kit"
        ]

        columns = []
        for cc in COVERAGE_COLUMNS:
            if getattr(self, cc) is not None:
                columns.append(cc)

        if not columns:
            msg = f"Couldn't find non-null coverage columns in QCExecSummary: {self.pk!r}"
            raise ValueError(msg)

        return columns

    @staticmethod
    def load_for_qc(_seqauto_run, qc, **kwargs):
        exec_summary_filename = QC.get_tsv_path_from_vcf(qc.vcf_file)
        exec_summary_data = load_exec_summary(QCExecSummary, exec_summary_filename)

        # Sanity check sample names match
        # TODO: Better name or way of storing this info than "aligned_pattern"?
        params = qc.get_params()
        expected_sample_name = params['sample_name']
        sample = exec_summary_data["sample"]
        if sample != expected_sample_name:
            msg = f"Exec summary file '{qc.path}' had sample name of '{sample}' while we expected '{expected_sample_name}'"
            raise ValueError(msg)

        exec_data = exec_summary_data["exec_data"]
        exec_summary = QCExecSummary.objects.create(qc=qc, sequencing_run=qc.sequencing_run,
                                                    data_state=DataState.COMPLETE, **exec_data)

        reference_range = exec_summary_data["reference_range"]
        if reference_range:
            reference_range["exec_summary"] = exec_summary
            ExecSummaryReferenceRange.objects.create(**reference_range)

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


class QCGeneCoverage(SeqAutoRecord):
    qc = models.OneToOneField(QC, on_delete=CASCADE)
    gene_coverage_collection = models.OneToOneField(GeneCoverageCollection, null=True, on_delete=CASCADE)

    @staticmethod
    def get_path_from_qc(qc):
        pattern = os.path.join(settings.SEQAUTO_QC_DIR_PATTERN, settings.SEQAUTO_QC_GENE_COVERAGE_PATTERN)
        return pattern % qc.get_params()

    def load_from_file(self, seqauto_run, **kwargs):
        enrichment_kit = None
        if self.qc:  # SeqAuto
            enrichment_kit = self.qc.sequencing_sample.enrichment_kit

        gene_coverage_collection = GeneCoverageCollection.objects.create(path=self.path,
                                                                         genome_build=self.qc.genome_build)
        self.gene_coverage_collection = gene_coverage_collection
        self.data_state = DataState.RUNNING
        self.save()

        try:
            warnings = gene_coverage_collection.load_from_file(enrichment_kit, **kwargs)
            if seqauto_run:
                for w in warnings:
                    SeqAutoMessage.objects.update_or_create(record=self, message=w, severity=LogLevel.WARNING,
                                                            defaults={"seqauto_run": seqauto_run})
            self.data_state = DataState.COMPLETE
        except FileNotFoundError:
            self.data_state = DataState.DELETED

        self.save()


@receiver(pre_delete, sender=QCGeneCoverage)
def gene_coverage_collection_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    try:
        if gcc := instance.gene_coverage_collection:
            instance.gene_coverage_collection = None  # To stop recursive deleting
            instance.save()
            gcc.delete()
    except:
        # Might fail due to GoldGeneCoverageCollection protecting it
        pass


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
    transcript_version = models.ForeignKey(TranscriptVersion, null=True, on_delete=SET_NULL)
    original_gene_symbol = models.TextField()
    original_transcript = models.TextField()
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


class JobScript(SeqAutoRecord):
    # Attach the SampleSheet/bam/vcf/qc via record so delete will cascade through
    seqauto_run = models.ForeignKey(SeqAutoRun, on_delete=CASCADE)
    out_file = models.TextField(null=True)
    file_type = models.CharField(max_length=1, choices=SequencingFileType.choices)
    record = models.ForeignKey(SeqAutoRecord, on_delete=CASCADE, related_name="+")
    job_id = models.TextField(null=True)
    job_status = models.CharField(max_length=1, choices=JobScriptStatus.choices, default=JobScriptStatus.CREATED)
    return_code = models.IntegerField(null=True)

    def get_record(self) -> SeqAutoRecord:
        return SeqAutoRecord.objects.get_subclass(pk=self.record_id)

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
        return f"{self.file_type}: {record_pk!r}"


@receiver(post_delete, sender=JobScript)
def post_delete_job_script(sender, instance, **kwargs):  # pylint: disable=unused-argument
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
    gcg_qs = GeneCoverageCollection.objects.filter(
        qcgenecoverage__qc__bam_file__unaligned_reads__sequencing_sample__sample_sheet__sequencingruncurrentsamplesheet__isnull=False)

    existing_count = 0
    existing_max_pk = 0
    existing_num = 0
    cache_key = f"20x_gene_coverage_{gene_symbol}_cov_{min_coverage}"
    if cached_count_num_and_max_pk := cache.get(cache_key):
        cached_count, cached_num_coverage_collections, cached_max_pk = cached_count_num_and_max_pk
        # Check that num is what we expect, ie none have been deleted/changed etc
        actual_num = gcg_qs.filter(pk__lte=cached_max_pk).count()
        if cached_num_coverage_collections == actual_num:
            logging.info("get_20x_gene_coverage - using cached values %d / %d / %d", *cached_count_num_and_max_pk)
            existing_count = cached_count
            existing_num = cached_num_coverage_collections
            existing_max_pk = cached_max_pk
            gcg_qs = gcg_qs.filter(pk__gt=existing_max_pk)
        else:
            logging.warning("get_20x_gene_coverage cache out of date (num records don't match actual)")

    if new_num := gcg_qs.count():
        logging.info("get_20x_gene_coverage - retrieving counts on %d new collections", new_num)
        # Using gene_coverage_collection__in=gcg_qs didn't save any time (SQL must evaluate inner query last)
        qs = GeneCoverageCanonicalTranscript.objects.filter(gene_symbol=gene_symbol,
                                                            gene_coverage_collection__pk__gte=existing_max_pk,
                                                            percent_20x__gte=min_coverage)
        count = existing_count + qs.count()
        data = gcg_qs.aggregate(max_pk=Max("pk"))
        max_pk = data["max_pk"] or existing_max_pk
        num = existing_num + new_num
        logging.info("get_20x_gene_coverage setting cache values %d / %d / %d", count, num, max_pk)
        cache.set(cache_key, (count, num, max_pk), timeout=30 * DAY_SECS)
    else:
        count = existing_count
    return count


def get_variant_caller_from_vcf_file(vcf_path):
    variant_caller, version = get_variant_caller_and_version_from_vcf(vcf_path)
    if variant_caller is None:
        variant_caller = "Unknown Variant Caller"
        version = -1

    return VariantCaller.objects.get_or_create(name=variant_caller, version=version)[0]
