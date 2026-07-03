import logging
import os
import pathlib
import re
from datetime import datetime
from functools import cached_property
from typing import Optional

from cache_memoize import cache_memoize
from django.conf import settings
from django.contrib.auth.models import User, Group
from django.contrib.postgres.fields import DecimalRangeField
from django.core.cache import cache
from django.db import models
from django.db.models import Max
from django.db.models.deletion import SET_NULL, CASCADE, PROTECT
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.urls.base import reverse
from django.utils.timezone import make_aware
from django_extensions.db.models import TimeStampedModel

from genes.models import GeneListCategory, CustomTextGeneList, GeneList, GeneCoverageCollection, \
    Transcript, GeneSymbol, SampleGeneList, TranscriptVersion, GeneCoverageCanonicalTranscript, ActiveSampleGeneList
from library.constants import DAY_SECS
from library.genomics.vcf_utils import get_variant_caller_and_version_from_vcf
from library.preview_request import PreviewModelMixin
from library.utils import sorted_nicely
from library.utils.file_utils import name_from_filename
from patients.models import FakeData, Patient
from seqauto.illumina.illumina_sequencers import SEQUENCING_RUN_REGEX
from seqauto.models.models_enums import DataGeneration, SequencerRead, PairedEnd
from seqauto.models.models_sequencing import Sequencer, EnrichmentKit, Experiment
from seqauto.models.models_software import Aligner, VariantCaller
from seqauto.signals.signals_list import sequencing_run_sample_sheet_created_signal
from snpdb.models import VCF, Sample, GenomeBuild, DataState, InheritanceManager, Wiki
from snpdb.models.models_enums import ImportStatus


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
        self.is_valid = True
        super().save(*args, **kwargs)

    def validate(self) -> bool:
        return True

    @property
    def last_modified_datetime(self):
        return make_aware(datetime.fromtimestamp(self.file_last_modified))

    @staticmethod
    def get_file_last_modified(filename):
        path = pathlib.Path(filename)
        return path.stat().st_mtime


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

    def get_current_sample_sheet(self):
        try:
            return self.sequencingruncurrentsamplesheet.sample_sheet
        except Exception:
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

    def get_flowcell_id(self) -> str:
        run_name = self.get_original_illumina_sequencing_run(self.name)
        return run_name.split("_")[-1]

    @staticmethod
    def get_date_from_name(name) -> Optional[datetime.date]:
        date = None
        if m := re.match(r".*_?([12]\d{5})_", name):
            date_str = m.group(1)
            dt = datetime.strptime(date_str, "%y%m%d")
            date = make_aware(dt).date()
        return date

    @property
    def original_sequencing_run(self):
        return self.get_original_illumina_sequencing_run(self.name)

    @cache_memoize(DAY_SECS, args_rewrite=lambda p: (p.pk, ))
    def get_params(self):
        """ This allows chaining down names etc - in case a level changes it, will cascade down """
        params = {
            "sequencing_run": self.name,
            "sequencing_run_dir": self.path,
            "flowcell_id": self.get_flowcell_id(),
            "original_sequencing_run": self.original_sequencing_run,
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

    def get_old_sample_sheets(self):
        current_sample_sheet = self.get_current_sample_sheet()
        return SampleSheet.objects.filter(sequencing_run=self).exclude(pk=current_sample_sheet.pk)

    @property
    def is_data_out_of_date_from_current_sample_sheet(self):
        try:
            current_sample_sheet = self.get_current_sample_sheet()
        except SequencingRunCurrentSampleSheet.DoesNotExist:
            return False

        for joint_called_vcf in current_sample_sheet.jointcalledvcf_set.all():
            if joint_called_vcf.needs_to_be_linked():
                return True

        old_sample_sheets = self.get_old_sample_sheets()
        illuminate_qc = IlluminaFlowcellQC.objects.filter(sample_sheet=current_sample_sheet)
        old_illuminate_qc = IlluminaFlowcellQC.objects.filter(sample_sheet__in=old_sample_sheets)
        old_joint_called_vcfs = JointCalledVCF.objects.filter(sample_sheet__in=old_sample_sheets)
        old_sample_links = SampleFromSequencingSample.objects.filter(sequencing_sample__sample_sheet__in=old_sample_sheets)
        old_unaligned_reads = UnalignedReads.get_for_old_sample_sheets(self)
        return any([not illuminate_qc.exists() and old_illuminate_qc.exists(),
                    old_joint_called_vcfs.exists(),
                    old_unaligned_reads.exists(),
                    old_sample_links.exists()])

    @staticmethod
    def get_external_links_for(name: str, date: Optional[datetime.date],
                               enrichment_kit_name: Optional[str]) -> list[tuple[str, str]]:
        """ Returns (label, url) tuples for external systems configured in
            SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS that apply to a run with these fields.
            Used from both the detail page (via get_external_links) and the runs grid. """
        links = []
        params = None
        for cfg in settings.SEQAUTO_SEQUENCING_RUN_EXTERNAL_LINKS:
            min_date_str = cfg.get("min_date")
            if min_date_str:
                min_date = datetime.strptime(min_date_str, "%Y-%m-%d").date()
                if not date or date < min_date:
                    continue
            excluded = cfg.get("exclude_enrichment_kits") or []
            if excluded and enrichment_kit_name:
                kit_name_lower = enrichment_kit_name.lower()
                if any(ex.lower() == kit_name_lower for ex in excluded):
                    continue
            if params is None:
                original = SequencingRun.get_original_illumina_sequencing_run(name)
                params = {
                    "sequencing_run": name,
                    "original_sequencing_run": original,
                    "flowcell_id": original.split("_")[-1],
                    "enrichment_kit": enrichment_kit_name or "",
                }
            links.append((cfg["label"], cfg["url_pattern"] % params))
        return links

    def get_external_links(self) -> list[tuple[str, str]]:
        kit_name = self.enrichment_kit.name if self.enrichment_kit else None
        return self.get_external_links_for(self.name, self.date, kit_name)

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

    def get_version_string(self):
        date = "TODO"
        return f"{self.hash}/{date}"

    def set_as_current_sample_sheet(self, sequencing_run, created=False):
        if created:
            sequencing_run_sample_sheet_created_signal.send(sender=os.path.basename(__file__),
                                                            sample_sheet=self)

        # Make sure SequencingRunCurrentSampleSheet is set to what the API sent us
        try:
            # Update existing
            current_ss = sequencing_run.sequencingruncurrentsamplesheet
            on_disk_not_current = current_ss.sample_sheet != self
            if on_disk_not_current:
                from seqauto.sequencing_files.sample_sheet import current_sample_sheet_changed
                current_sample_sheet_changed(current_ss, self)

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
    barcode = models.TextField()  # historically we stored 'index' in here. Now we store index1|index2
    enrichment_kit = models.ForeignKey(EnrichmentKit, null=True, on_delete=CASCADE)
    is_control = models.BooleanField(default=False)
    failed = models.BooleanField(default=False)
    automatically_process = models.BooleanField(default=True)

    @cache_memoize(DAY_SECS, args_rewrite=lambda p: (p.pk, ))
    def get_params(self):
        params = self.sample_sheet.get_params()
        indexes = self.barcode.split("|")
        index = indexes[0]
        if len(indexes) > 1:
            index2 = indexes[1]
        else:
            index2 = "NO_INDEX2_STORED"
        lane_num = self.lane or 1
        params.update({"lane": lane_num,
                       "lane_code": f"L{lane_num:03d}",
                       "sample_number": self.sample_number,
                       "sample_id": self.sample_id,
                       "sample_name": self.sample_name,
                       "sample_project": self.sample_project or '',
                       "sample_name_underscores": self.sample_name.replace("-", "_"),
                       "barcode": index,
                       "index": index,
                       "index2": index2})
        if len(indexes) == 2:
            params["index2"] = indexes[1]

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
            except Exception:
                logging.error("Wasn't exactly 1 bam_file for unaligned_reads %s", unaligned_reads)
        except Exception:
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

    def get_variant_caller(self):
        return self.variant_caller


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

    def get_params(self):
        params = self.sequencing_sample.get_params()
        params['fastq'] = self.path
        return params

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

    @cache_memoize(DAY_SECS, args_rewrite=lambda p: (p.pk, ))
    def get_params(self):
        fastq_params = self.sequencing_sample.get_params()
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
        if self.fastq_r2 and self.fastq_r1.data_state != self.fastq_r2.data_state:
            data_state = f"R1: {self.fastq_r1.get_data_state_display()}, R2: {self.fastq_r2.get_data_state_display()}"
        else:
            data_state = self.fastq_r1.get_data_state_display()
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

    def __str__(self):
        return f"Flagstats ({self.get_data_state_display()}) for {self.bam_file}"


def get_seqauto_user():
    user, created = User.objects.get_or_create(username=settings.SEQAUTO_USER)
    if created:
        if settings.SEQAUTO_GROUP:
            seqauto_group = Group.objects.get_or_create(name=settings.SEQAUTO_GROUP)
            user.groups.add(seqauto_group)
    return user


class SingleSampleVCF(SeqAutoRecord):
    """ Single-sample VCFs from the file system, tied to exactly one BamFile """
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
        except Exception:
            up = None
        return up

    @property
    def sample_sheet(self) -> SampleSheet:
        return self.bam_file.unaligned_reads.sequencing_sample.sample_sheet

    def __str__(self):
        return f"VCF {self.name} ({self.get_data_state_display()})"


class JointCalledVCF(SeqAutoRecord):
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
        except Exception:
            up = None
        return up

    @cached_property
    def vcf(self) -> Optional[VCF]:
        try:
            return VCF.objects.get(uploadedvcf__uploaded_file__path=self.path)
        except VCF.DoesNotExist:
            return None

    @property
    def sample_count(self) -> Optional[int]:
        if vcf := self.vcf:
            return vcf.sample_set.count()
        return None

    @property
    def is_full_sheet(self) -> Optional[bool]:
        count = self.sample_count
        if count is None:
            return None
        return count == self.sample_sheet.sequencingsample_set.count()

    def needs_to_be_linked(self):
        try:
            _ = self.backendvcf  # if no exception we can link it
            try:
                _ = self.backendvcf.vcf.vcffromsequencingrun
            except VCFFromSequencingRun.DoesNotExist:
                return True
        except Exception:
            pass
        return False

    def __str__(self):
        num_samples = self.sample_sheet.sequencingsample_set.count()
        return f"{self.variant_caller} JointCalledVCF ({self.pk}) for {self.sequencing_run} ({num_samples} samples)"


class QC(SeqAutoRecord):
    """ This holds together all of the other QC objects

        Historically, we regarded "QC" as the exec stats file, then we needed to load a whole bunch of QC so those
        became new models that hung off "QC" which is just used to group things

        So because of this, QC is dependent on "path". Probably the "right" thing to do is to change QC to not be a
        SeqAuto object, as it's not actually a real object, just a way to link models.

        TODO: Make this not a SeqAutoRecord. Perhaps make this unique_together w/bam+vcf?
    """
    bam_file = models.ForeignKey(BamFile, on_delete=CASCADE)
    vcf_file = models.ForeignKey(SingleSampleVCF, on_delete=CASCADE)

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

    def link_samples_if_exist(self, force_active=False):
        # SampleFromSequencingSample is only done after VCF import, so this may not be linked yet.
        for sfss in self.sequencing_sample.samplefromsequencingsample_set.all():
            sample = sfss.sample
            self.create_and_assign_sample_gene_list(sample, force_active=force_active)  # Also saves

    def create_and_assign_sample_gene_list(self, sample: Sample, force_active=False):
        logging.info("QCGeneList: %d - create_and_assign_sample_gene_list for %s", self.pk, sample)
        # On create signal, 'sample_gene_list_created' sets to active gene list
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
    # gene_coverage_collection is populated in load_from_file, set_null as needs to be deleted/reset when reloaded
    gene_coverage_collection = models.OneToOneField(GeneCoverageCollection, null=True, on_delete=SET_NULL)

    @staticmethod
    def get_path_from_qc(qc):
        pattern = os.path.join(settings.SEQAUTO_QC_DIR_PATTERN, settings.SEQAUTO_QC_GENE_COVERAGE_PATTERN)
        return pattern % qc.get_params()


@receiver(pre_delete, sender=QCGeneCoverage)
def gene_coverage_collection_pre_delete_handler(sender, instance, **kwargs):  # pylint: disable=unused-argument
    try:
        if gcc := instance.gene_coverage_collection:
            instance.gene_coverage_collection = None  # To stop recursive deleting
            instance.save()
            gcc.delete()
    except Exception:
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
            name += f" error: {self.error_exception}"
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
        except Exception:
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


def get_samples_by_sequencing_sample(sample_sheet, vcf):
    sequencing_samples_by_name = sample_sheet.get_sequencing_samples_by_name()

    def clean_sample_name(s):
        return s.upper().replace("-", "_")

    samples = list(vcf.sample_set.all())
    potential_samples_by_name = {
        s.name: s for s in samples
    }
    # If there is a single sample, use VCF name (work around SpliceGirls samples called "SAMPLE")
    if len(samples) == 1:
        potential_samples_by_name[vcf.name] = samples[0]

    samples_by_sequencing_sample = {}
    for sample_name, sample in potential_samples_by_name.items():
        # Do a startswith match rather than hash lookup as it's less strict (diff naming conventions etc)
        for sequencing_sample_name, sequencing_sample in sequencing_samples_by_name.items():
            cleaned_sample_name = clean_sample_name(sample_name)
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
