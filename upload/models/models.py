import logging
import os
import re
import shutil
from collections import namedtuple
from typing import Dict, List, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ObjectDoesNotExist, PermissionDenied
from django.db import models, transaction
from django.db.models.aggregates import Max
from django.db.models.deletion import CASCADE, SET_NULL
from django.db.models.query import QuerySet
from django.db.models.signals import pre_delete
from django.dispatch.dispatcher import receiver
from django.utils import timezone
from django_extensions.db.models import TimeStampedModel
from model_utils.managers import InheritanceManager

from eventlog.models import create_event
from library.django_utils.django_file_system_storage import PrivateUploadStorage
from library.django_utils.django_file_utils import get_import_processing_dir
from library.enums.log_level import LogLevel
from library.file_utils import file_md5sum, mk_path
from library.log_utils import get_traceback, report_message
from seqauto.models import VCFFile, SampleSheetCombinedVCFFile, get_samples_by_sequencing_sample, VariantCaller
from snpdb.import_status import set_vcf_and_samples_import_status
from snpdb.models import VCF, Variant, SoftwareVersion, GenomeBuild
from snpdb.models.models_enums import ImportSource, ProcessingStatus, ImportStatus
from snpdb.tasks.soft_delete_tasks import soft_delete_vcfs
from upload.models.models_enums import UploadedFileTypes, VCFPipelineStage, \
    UploadStepTaskType, TimeFilterMethod, VCFImportInfoSeverity, UploadStepOrigin
from variantgrid.celery import app


class UploadedFile(TimeStampedModel):
    path = models.TextField(null=True)
    uploaded_file = models.FileField(storage=PrivateUploadStorage(),
                                     max_length=256, null=True)
    md5_hash = models.CharField(max_length=32, null=True)
    file_type = models.CharField(max_length=1, choices=UploadedFileTypes.choices, null=True)
    import_source = models.CharField(max_length=1, choices=ImportSource.choices)
    name = models.TextField()
    user = models.ForeignKey(User, on_delete=CASCADE)
    visible = models.BooleanField(default=True)  # Set False to hide on upload page

    @property
    def exists(self):
        if self.uploaded_file:
            return os.path.exists(self.get_filename())
        return False

    @property
    def size(self) -> int:
        if self.import_source == ImportSource.WEB_UPLOAD:
            return self.uploaded_file.size
        return os.stat(self.get_filename()).st_size

    def check_can_view(self, user):
        if not (user.is_superuser or self.user == user):
            msg = f"You do not have permission to access UploadedFile pk={self.pk}"
            raise PermissionDenied(msg)

    def get_file(self):
        return open(self.get_filename(), "rb")

    def get_filename(self):
        if self.import_source == ImportSource.WEB_UPLOAD:
            filename = self.uploaded_file.path
        else:
            filename = self.path
        return filename

    def store_md5_hash(self):
        self.md5_hash = file_md5sum(self.get_filename())
        self.save()

    def __str__(self):
        description = f"{self.name} ({self.created})"
        return description


class PipelineFailedJobTerminateEarlyException(Exception):
    """ Throw this if pipeline fails elsewhere, and want to exit early but not log exception """
    pass


class SkipUploadStepException(Exception):
    """ Throw this to set task to skipped """
    pass


class UploadPipeline(models.Model):
    INITIAL_PROGRESS_STATUS = "Processing"
    status = models.CharField(max_length=1, choices=ProcessingStatus.choices, default=ProcessingStatus.CREATED)
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    items_processed = models.BigIntegerField(null=True)
    processing_seconds_wall_time = models.IntegerField(null=True)
    processing_seconds_cpu_time = models.IntegerField(null=True)
    progress_status = models.TextField(null=True)
    progress_percent = models.FloatField(default=0)
    celery_task = models.CharField(max_length=36, null=True)

    @property
    def file_type(self):
        return self.uploaded_file.file_type

    @property
    def get_file_type_display(self):
        return self.uploaded_file.get_file_type_display()

    def get_max_step_sort_order(self):
        qs = self.uploadstep_set.exclude(pipeline_stage_dependency=VCFPipelineStage.FINISH)
        data = qs.aggregate(Max("sort_order"))
        return data["sort_order__max"] or 0

    def get_pipeline_processing_dir(self):
        return get_import_processing_dir(self.pk)

    def get_pipeline_processing_subdir(self, subdir):
        pipeline_processing_dir = self.get_pipeline_processing_dir()
        sd = os.path.join(pipeline_processing_dir, subdir)
        mk_path(sd)
        return sd

    def remove_processing_files(self):
        pipeline_processing_dir = self.get_pipeline_processing_dir()
        logging.info("*** Deleting files for pipeline %d - '%s'", self.pk, pipeline_processing_dir)
        shutil.rmtree(pipeline_processing_dir)

    def start(self):
        logging.debug("upload_pipeline.start()")
        self.status = ProcessingStatus.PROCESSING
        self.items_processed = 0
        self.progress_status = UploadPipeline.INITIAL_PROGRESS_STATUS
        self.progress_percent = 0
        self.save()

        user = self.uploaded_file.user
        create_event(user, f"import_{self.file_type}_start")

    def success(self, items_processed=None, processing_seconds_wall_time=None, processing_seconds_cpu_time=None):
        if self.progress_status == UploadPipeline.INITIAL_PROGRESS_STATUS:
            self.progress_status = 'Success'
        self.status = ProcessingStatus.SUCCESS
        self.items_processed = items_processed
        self.processing_seconds_wall_time = processing_seconds_wall_time
        self.processing_seconds_cpu_time = processing_seconds_cpu_time
        self.progress_percent = 100.0
        self.save()

        user = self.uploaded_file.user
        create_event(user, f"import_{self.file_type}_success")
        if settings.VCF_IMPORT_DELETE_TEMP_FILES_ON_SUCCESS:
            self.remove_processing_files()

    def error(self, error_message):
        logging.error("upload_pipeline.error(%s, %s)", self, error_message)

        self.status = ProcessingStatus.ERROR
        self.progress_status = f"Error: {error_message}"
        self.save()

        if vcf := self.vcf:
            set_vcf_and_samples_import_status(vcf, ImportStatus.ERROR)

        user = self.uploaded_file.user
        name = f"import_{self.file_type}_failed"
        filename = self.uploaded_file.get_filename()
        create_event(user, name, error_message, filename=filename, severity=LogLevel.ERROR)
        message = f"UploadPipeline {self.pk} failed. " \
            f"Filename: {self.get_file_type_display} Error: {error_message}"
        report_message(message, level='error')

    @property
    def vcf(self):
        if self.uploaded_file.file_type == UploadedFileTypes.VCF:
            try:
                return self.uploaded_file.uploadedvcf.vcf
            except (UploadedVCF.DoesNotExist, VCF.DoesNotExist) as _:
                pass
        return None

    def get_errors(self, hide_accepted=True) -> List:
        errors = []
        if self.status == 'E':
            # Make a fake one for template
            VCFImportError = namedtuple('VCFImportError', ['message', 'has_more_details'])
            vii = VCFImportError(message=f"This file failed to import due to: {self.progress_status}.",
                                 has_more_details=False)
            errors.append(vii)
        errors.extend(self._get_vcf_import_info(VCFImportInfoSeverity.ERROR, hide_accepted=hide_accepted))
        return errors

    def get_warnings(self, hide_accepted=True, include_vcf=True) -> List:
        warnings = self._get_vcf_import_info(VCFImportInfoSeverity.WARNING, hide_accepted=hide_accepted)
        if include_vcf:
            if vcf := self.vcf:
                warnings.extend(vcf.get_warnings())
        return warnings

    def _get_vcf_import_info(self, severity, hide_accepted=True):
        kwargs = {"upload_step__upload_pipeline": self,
                  "severity": severity}
        if hide_accepted:
            kwargs["accepted_date__isnull"] = True
        vcf_import_info = []
        for import_info in VCFImportInfo.objects.filter(**kwargs).select_subclasses():
            msg = import_info.message
            if msg:
                vcf_import_info.append(import_info)
        return vcf_import_info

    @property
    def genome_build(self):
        """ returns GenomeBuild based on uploaded file type
            throws ValueError if can't retrieve it """

        # TODO: bit of a hack - do this more OO?
        uploaded_file = self.uploaded_file
        file_type = uploaded_file.file_type

        genome_build = None
        if file_type == UploadedFileTypes.VCF:
            genome_build = uploaded_file.uploadedvcf.vcf.genome_build
        elif file_type == UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY:
            # This can be either a manually entered variants or from variant classification
            try:
                genome_build = uploaded_file.uploadedmanualvariantentrycollection.collection.genome_build
            except ObjectDoesNotExist:
                try:
                    genome_build = uploaded_file.uploadedclassificationimport.classification_import.genome_build
                except ObjectDoesNotExist:
                    pass
        elif file_type == UploadedFileTypes.LIFTOVER:
            genome_build = uploaded_file.uploadedliftover.liftover.genome_build
        elif file_type == UploadedFileTypes.CLINVAR:
            genome_build = uploaded_file.uploadedclinvarversion.clinvar_version.genome_build
        elif file_type == UploadedFileTypes.VARIANT_TAGS:
            genome_build = uploaded_file.uploadedvarianttags.variant_tags_import.genome_build

        if genome_build is None:
            msg = f"Don't know how to get GenomeBuild for UploadedFile type '{uploaded_file.get_file_type_display()}'"
            raise ValueError(msg)

        return genome_build

    def __str__(self):
        return f"UploadPipeline {self.pk} ({self.file_type}: {self.status})"


class ToolVersion(SoftwareVersion):
    pass


class UploadStep(models.Model):
    """ This is only used for VCF multi-stage imports """
    CREATE_UNKNOWN_LOCI_AND_VARIANTS_TASK_NAME = "Create Unknown Loci and Variants"
    PREPROCESS_VCF_NAME = "Preprocess VCF"
    PROCESS_VCF_TASK_NAME = "Process VCF File"

    name = models.TextField()
    upload_pipeline = models.ForeignKey(UploadPipeline, on_delete=CASCADE)
    parent_upload_step = models.ForeignKey('self', null=True, related_name='substep_set', on_delete=CASCADE)
    sort_order = models.IntegerField()
    status = models.CharField(max_length=1, choices=ProcessingStatus.choices, default=ProcessingStatus.CREATED)
    # origin is used to indicate whether a step was from a import task factory where it will
    # be re-generated upon retry-import, or manually created (will be left alone and re-executed)
    origin = models.CharField(max_length=1, choices=UploadStepOrigin.choices, default=UploadStepOrigin.IMPORT_TASK_FACTORY)
    items_to_process = models.IntegerField(null=True)
    items_processed = models.IntegerField(null=True)
    error_message = models.TextField()
    input_upload_step = models.ForeignKey("self", null=True, on_delete=CASCADE)
    input_filename = models.TextField()
    output_filename = models.TextField()
    start_date = models.DateTimeField(null=True)
    end_date = models.DateTimeField(null=True)
    task_type = models.CharField(max_length=1, choices=UploadStepTaskType.choices)
    pipeline_stage = models.CharField(max_length=1, choices=VCFPipelineStage.choices, null=True)
    pipeline_stage_dependency = models.CharField(max_length=1, choices=VCFPipelineStage.choices, null=True)
    script = models.TextField()
    child_script = models.TextField(null=True)
    import_variant_table = models.TextField(null=True, blank=True)
    celery_task = models.CharField(max_length=36, null=True)
    tool_version = models.ForeignKey(ToolVersion, null=True, on_delete=CASCADE)
    output_text = models.TextField(null=True)

    @property
    def uploaded_file(self) -> UploadedFile:
        return self.upload_pipeline.uploaded_file

    @property
    def processing_seconds(self):
        if self.start_date and self.end_date:
            seconds = (self.end_date - self.start_date).total_seconds()
        else:
            seconds = None
        return seconds

    @property
    def genome_build(self):
        return self.upload_pipeline.genome_build

    def start(self):
        if self.status != ProcessingStatus.CREATED:
            msg = f"{self.pk} called start() with process {self.get_status_display()} should be CREATED"
            raise ValueError(msg)

        self.start_date = timezone.now()
        self.status = ProcessingStatus.PROCESSING
        self.items_processed = 0

    def error_exception(self, e):
        status = ProcessingStatus.ERROR
        self.status = status
        self.error_message = str(e) + "\n" + get_traceback()

        qs = UploadPipeline.objects.filter(pk=self.upload_pipeline.pk)
        progress_status = f"Failure in import step id: {self.pk}, ({self.error_message})"
        qs.update(status=status, progress_status=progress_status)

    @transaction.atomic()
    def mark_timed_out(self, user: Optional[User] = None):
        self.status = ProcessingStatus.TIMED_OUT
        if user:
            self.error_message = f"Marked as timed out by {user.username}"
        else:
            self.error_message = "Timed out"

        self.save()
        qs = UploadPipeline.objects.filter(pk=self.upload_pipeline.pk)
        progress_status = f"Failure in import step id: {self.pk}, ({self.error_message})"
        qs.update(status=self.status, progress_status=progress_status)

    def get_uploaded_vcf(self) -> 'UploadedVCF':
        return self.uploaded_file.uploadedvcf

    def get_multi_input_files_and_records(self):
        input_filenames_and_records = []
        if self.input_upload_step:
            for upsmfo in self.input_upload_step.uploadstepmultifileoutput_set.all():
                input_filenames_and_records.append((upsmfo.output_filename, upsmfo.items_to_process))

        return input_filenames_and_records

    def pipeline_inserted_unknown_variants(self):
        steps_qs = self.upload_pipeline.uploadstep_set.all()
        create_variants_step = steps_qs.filter(name=UploadStep.CREATE_UNKNOWN_LOCI_AND_VARIANTS_TASK_NAME)
        return create_variants_step.exists()

    def launch_task(self, task_clazz):
        logging.info("task_clazz: %s", task_clazz)
        task = task_clazz.si(self.pk, 0)
        result = task.apply_async()

        celery_task = result.id
        UploadStep.objects.filter(pk=self.pk).update(celery_task=celery_task)

    def close_sub_steps(self):
        for sub_step in self.substep_set.all():
            sub_step.status = self.status
            sub_step.end_date = self.end_date
            sub_step.save()

    def __str__(self):
        return f"UploadStep: {self.name} (pk={self.pk}) ({self.upload_pipeline.pk}/{self.status})"


class UploadStepMultiFileOutput(models.Model):
    """ A step may split a file up into multiple files (eg per chrom)  """
    upload_step = models.ForeignKey(UploadStep, on_delete=CASCADE)
    output_filename = models.TextField()
    items_to_process = models.IntegerField(null=True)


class VCFImporter(models.Model):
    """ Keep track of how the processing is done, in case we need to audit/fix import bugs later """
    name = models.TextField()
    version = models.IntegerField()
    vcf_parser = models.TextField()
    vcf_parser_version = models.TextField()
    code_git_hash = models.TextField()

    class Meta:
        unique_together = ("name", "version", "vcf_parser", "vcf_parser_version", "code_git_hash")

    def __str__(self):
        return f"{self.name} (v.{self.version}). Git: {self.code_git_hash}. Uses {self.vcf_parser} (v.{self.vcf_parser_version})"


# Note: Probably the best thing to do is to make the below subclasses UploadedFile but migration might be a hassle
class UploadedVCF(models.Model):
    uploaded_file = models.OneToOneField(UploadedFile, on_delete=CASCADE)
    vcf = models.OneToOneField(VCF, null=True, on_delete=CASCADE)
    upload_pipeline = models.OneToOneField(UploadPipeline, null=True, on_delete=SET_NULL)
    # Keep track of highest KNOWN variant, so we can quickly check whether it has been annotated.
    max_variant = models.ForeignKey(Variant, null=True, on_delete=SET_NULL)
    vcf_importer = models.ForeignKey(VCFImporter, null=True, on_delete=SET_NULL)

    def get_data(self) -> VCF:
        return self.vcf

    def get_upload_context(self) -> Dict:
        """ Dict for displaying JFU upload widget """
        context = {}
        if self.vcf:
            context["import_status"] = self.vcf.import_status
            samples = []
            for sample in self.vcf.sample_set.all():
                samples.append({"sample_id": sample.pk,
                                "name": sample.name})
            context["samples"] = samples

        return context

    def __str__(self):
        description = f"UploadedVCF for {self.uploaded_file}"
        if self.vcf:
            description += f" (proj: {self.vcf})"
        return description


@receiver(pre_delete, sender=UploadedVCF)
def pre_delete_uploaded_vcf(sender, instance, *args, **kwargs):
    if vcf := instance.vcf:
        #logging.debug("Deleting associated vcf for uploaded VCF")
        if vcf.import_status not in ImportStatus.DELETION_STATES:
            soft_delete_vcfs(vcf.user, vcf.pk)


class UploadedVCFPendingAnnotation(models.Model):
    uploaded_vcf = models.OneToOneField(UploadedVCF, on_delete=CASCADE)
    created = models.DateTimeField(auto_now_add=True)
    finished = models.DateTimeField(null=True)
    schedule_pipeline_stage_steps_celery_task = models.CharField(max_length=36, null=True)

    def attempt_schedule_annotation_stage_steps(self, lowest_unannotated_variant_id):
        if self.uploaded_vcf.max_variant is None or self.uploaded_vcf.max_variant_id < lowest_unannotated_variant_id:
            logging.info("UploadedVCF %s all variants are annotated", self.uploaded_vcf)
            self.finished = timezone.now()
            self.save()

            # Has to be a string to get around circular imports
            CELERY_TASK_NAME = 'upload.tasks.vcf.import_vcf_step_task.schedule_pipeline_stage_steps'
            logging.info("Executing: %s", CELERY_TASK_NAME)
            args = (self.uploaded_vcf.upload_pipeline.pk, VCFPipelineStage.ANNOTATION_COMPLETE)
            result = app.send_task(CELERY_TASK_NAME, args=args)
            self.schedule_pipeline_stage_steps_celery_task = result.id
            self.save()
        else:
            logging.info("UploadedVCF %s still waiting for annotation", self.uploaded_vcf)


class BackendVCF(models.Model):
    """ Link between UploadedVCF (upload) and Filesystem VCF (SeqAuto) """
    uploaded_vcf = models.OneToOneField(UploadedVCF, on_delete=CASCADE)
    vcf_file = models.OneToOneField(VCFFile, null=True, on_delete=CASCADE)
    combo_vcf = models.OneToOneField(SampleSheetCombinedVCFFile, null=True, on_delete=CASCADE)

    @property
    def filesystem_vcf(self):
        vcfs = [self.combo_vcf, self.vcf_file]
        if all(vcfs):
            msg = f"{self} has both vcf_file and combo_vcf set"
            raise ValueError(msg)
        if self.vcf_file:
            return self.vcf_file
        if self.combo_vcf:
            return self.combo_vcf
        raise ValueError(f"{self} has neither 'vcf_file' or 'combo_vcf' set")

    @property
    def sample_sheet(self):
        return self.filesystem_vcf.sample_sheet

    @property
    def user(self):
        return self.uploaded_vcf.uploaded_file.user

    @property
    def vcf(self):
        return self.uploaded_vcf.vcf

    @property
    def variant_caller(self) -> VariantCaller:
        record = self.vcf_file or self.combo_vcf
        return record.variant_caller

    def get_samples_by_sequencing_sample(self):
        return get_samples_by_sequencing_sample(self.filesystem_vcf.sample_sheet, self.vcf)

    def __str__(self):
        backend = self.vcf_file or self.combo_vcf or 'None'
        return f"BackendVCF: backend: {backend}, uploaded_vcf: {self.uploaded_vcf}"


class VCFImportInfo(models.Model):
    objects = InheritanceManager()
    severity = models.CharField(max_length=1, choices=VCFImportInfoSeverity.choices, default=VCFImportInfoSeverity.WARNING)
    upload_step = models.ForeignKey(UploadStep, on_delete=CASCADE)
    accepted_date = models.DateTimeField(null=True)
    has_more_details = False

    @property
    def message(self):
        raise NotImplementedError()

    def __str__(self):
        return f"{self.get_severity_display()}: {self.message}"


class SimpleVCFImportInfo(VCFImportInfo):
    ANNOTATION_SKIPPED = "annotation_skipped"

    type = models.TextField()
    has_more_details = models.BooleanField(default=False)
    message_string = models.TextField()

    @property
    def message(self):
        return self.message_string


class ModifiedImportedVariants(VCFImportInfo):
    has_more_details = True

    @property
    def message(self):
        qs = self.modifiedimportedvariant_set.filter(old_variant__isnull=False)
        num_normalised = qs.count()
        if num_normalised:
            msg = f"{num_normalised} variants normalised during import."
        else:
            msg = None
        return msg


class ModifiedImportedVariant(models.Model):
    """ Keep track of variants that were modified during import pre-processing by vt,
        so people can find out why a variant they expected didn't turn up.  """
    VT_OLD_VARIANT_PATTERN = re.compile(r"^([^:]+):(\d+):([^/]+)/([^/]+)$")

    import_info = models.ForeignKey(ModifiedImportedVariants, on_delete=CASCADE, null=True)
    variant = models.ForeignKey(Variant, on_delete=CASCADE)
    # OLD_MULTIALLELIC from vt: @see https://genome.sph.umich.edu/wiki/Vt#Decompose
    old_multiallelic = models.TextField(null=True)
    # OLD_VARIANT from vt: @see https://genome.sph.umich.edu/wiki/Vt#Normalization
    old_variant = models.TextField(null=True)
    old_variant_formatted = models.TextField(null=True)  # consistently format for retrieval

    @property
    def genome_build(self):
        return self.import_info.upload_step.genome_build

    @staticmethod
    def format_old_variant(old_variant: str, genome_build: GenomeBuild) -> List[str]:
        """ We need consistent formatting (case and use of chrom) so we can retrieve it easily.
            May return multiple values """
        formatted_old_variants = []
        for ov in ModifiedImportedVariant._split_old_variant(old_variant):
            variant_tuple = Variant.get_tuple_from_string(ov, genome_build,
                                                          regex_pattern=ModifiedImportedVariant.VT_OLD_VARIANT_PATTERN)
            formatted_old_variants.append(ModifiedImportedVariant.get_old_variant_from_tuple(*variant_tuple))
        return formatted_old_variants

    @staticmethod
    def _split_old_variant(old_variant) -> List[str]:
        """ VT decompose writes OLD_VARIANT as comma separated, but if not decomposed (eg someone uploads already
            normalized) could be multi-alt that looks like 5:132240059:CT/CTT/T """
        old_variants = []
        for ov in old_variant.split(","):
            locus, *alts = ov.split("/")
            num_alts = len(alts)
            if num_alts == 1:
                # Quick most common path - no need to pull apart and re-assemble
                old_variants.append(ov)
            elif num_alts >= 2:
                for alt in alts:
                    old_variants.append(f"{locus}/{alt}")
            else:
                raise ValueError(f"Badly formatted INFO - 'OLD_VARIANT': '{ov}'")
        return old_variants

    @staticmethod
    def get_old_variant_from_tuple(chrom, position, ref, alt) -> str:
        return f"{chrom}:{int(position)}:{ref}/{alt}"

    @classmethod
    def get_upload_pipeline_unnormalized_variant(cls, upload_pipeline, chrom, position, ref, alt):
        """ throws DoesNotExist """
        old_variant = cls.get_old_variant_from_tuple(chrom, position, ref, alt)
        return cls.objects.get(import_info__upload_step__upload_pipeline=upload_pipeline,
                               old_variant_formatted=old_variant)

    @classmethod
    def get_variant_for_unnormalized_variant(cls, upload_pipeline, chrom, position, ref, alt) -> Variant:
        miv = cls.get_upload_pipeline_unnormalized_variant(upload_pipeline, chrom, position, ref, alt)
        return miv.variant

    @classmethod
    def get_variants_for_unnormalized_variant(cls, chrom, position, ref, alt) -> QuerySet:
        old_variant = cls.get_old_variant_from_tuple(chrom, position, ref, alt)
        return Variant.objects.filter(modifiedimportedvariant__old_variant_formatted=old_variant).distinct()

    @classmethod
    def get_variants_for_unnormalized_variant_any_alt(cls, chrom, position, ref) -> QuerySet:
        old_variant = cls.get_old_variant_from_tuple(chrom, position, ref, "")
        return Variant.objects.filter(modifiedimportedvariant__old_variant_formatted__startswith=old_variant).distinct()


class VCFSkippedContigs(VCFImportInfo):
    has_more_details = True

    @property
    def message(self):
        qs = self.vcfskippedcontig_set.all()
        d = dict(qs.values_list("contig", "num_skipped"))
        num_skipped_contigs = len(d)
        num_skipped_variants = sum(d.values())

        msg = None
        if num_skipped_contigs:
            msg = f"Import skipped {num_skipped_contigs} unknown contigs ({num_skipped_variants} total variants)."
        return msg


class VCFSkippedContig(models.Model):
    """ vcf_format_chrom filters unknown contigs and creates this """

    import_info = models.ForeignKey(VCFSkippedContigs, on_delete=CASCADE, null=True)
    contig = models.TextField()
    num_skipped = models.IntegerField()


class UploadSettings(models.Model):
    user = models.OneToOneField(User, on_delete=CASCADE)
    time_filter_method = models.CharField(max_length=1, choices=TimeFilterMethod.choices, default=TimeFilterMethod.RECORDS)
    time_filter_value = models.IntegerField(default=5)
    show_all = models.BooleanField(default=False)  # By default hide visible=False (liftover etc)

    def __str__(self):
        return f"User: {self.user}, {self.time_filter_value} {self.time_filter_method}"
