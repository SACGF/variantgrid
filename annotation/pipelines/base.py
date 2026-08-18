import abc
import logging
import os
from typing import Optional

from django.conf import settings
from django.db.models import QuerySet
from django.utils import timezone

from annotation.annotation_run_files import write_qs_to_vcf
from annotation.annotation_version_querysets import get_variants_qs_for_annotation
from annotation.models.models_enums import VariantAnnotationPipelineType
from library.utils.file_utils import mk_path_for_file
from snpdb.models import Variant


class AnnotationPipelineRunner(abc.ABC):
    """ One annotation pipeline's tool: which variants it takes, how it runs, and how its output is
        imported.

        The two abstract stages map to the two dispatcher lanes (#1649). annotate() executes on
        annotation_workers and leaves the run at ANNOTATION_COMPLETED with its output on disk;
        import_results() executes on db_workers and takes it to FINISHED. Everything around them - the
        task_id lock, lease renewal, reclaim, per-task paths, the attempt cap - belongs to the tasks in
        annotation.tasks.annotate_variants, so a new tool inherits all of it. """

    pipeline_type: VariantAnnotationPipelineType
    # VEP pipelines only want variants with no VariantAnnotation row yet - that is what "needs annotating"
    # means for them. Supplementary pipelines set this False and take every variant of their class in
    # range: they UPDATE rows a VEP pipeline wrote, and depends_on is what guarantees those rows exist.
    selects_unannotated: bool = True

    def get_variants_qs(self, annotation_run) -> QuerySet[Variant]:
        """ Variants in this run's range that this pipeline is responsible for. """
        range_lock = annotation_run.annotation_range_lock
        annotation_version = range_lock.version.get_any_annotation_version()
        return get_variants_qs_for_annotation(annotation_version,
                                              pipeline_type=self.pipeline_type,
                                              min_variant_id=range_lock.min_variant_id,
                                              max_variant_id=range_lock.max_variant_id,
                                              annotated=not self.selects_unannotated)

    def dump(self, annotation_run, dump_dir=None, task_token=None) -> int:
        """ Write this run's variants to a VCF and set the dump_* fields; returns the dump count.

            dump_dir lets the annotation_external --dump command (#1568) write elsewhere and stop before
            the tool runs. task_token (#1658) makes the local pipeline's dump path per-task so a reclaimed
            run can't collide with its zombie predecessor. """
        vcf_dump_filename = annotation_run.get_dump_filename(dump_dir=dump_dir, task_token=task_token)
        annotation_run.dump_start = timezone.now()
        annotation_run.vcf_dump_filename = vcf_dump_filename
        annotation_run.save()

        logging.info("Dumping %s variants for %s", self.pipeline_type, annotation_run)
        if os.path.exists(vcf_dump_filename):
            raise ValueError(f"Don't want to overwrite '{vcf_dump_filename}' which already exists!")
        mk_path_for_file(vcf_dump_filename)
        vcf_dump_count = write_qs_to_vcf(vcf_dump_filename, annotation_run.genome_build,
                                         self.get_variants_qs(annotation_run))

        annotation_run.dump_count = vcf_dump_count
        annotation_run.dump_end = timezone.now()
        annotation_run.save()
        return vcf_dump_count

    @abc.abstractmethod
    def annotate(self, annotation_run, lease_heartbeat=None):
        """ Dump this run's variants and run the tool over them, leaving the output on disk. """

    @abc.abstractmethod
    def import_results(self, annotation_run):
        """ Load the tool's output into the DB and set upload_end. """

    @abc.abstractmethod
    def get_output_paths(self, annotation_run, dump_filename) -> list[str]:
        """ Every file one attempt writes, derived from that attempt's dump stem.

            #1660: derivation rather than the persisted filename fields, which only reach the DB once the
            whole pipeline has finished - so a run reclaimed before then has files on disk its row cannot
            name. The dump stem is persisted before the tool starts, which makes it the one key that
            always names the rest. """

    @abc.abstractmethod
    def tool_finished(self, annotation_run, dump_filename) -> bool:
        """ Whether the tool ran to completion for this attempt AND its output is still on disk, judged
            from the files alone.

            Drives the reclaim's keep-or-scrub decision (_reset_run_for_redispatch): a run whose expensive
            output is complete resumes at the import step instead of re-running the tool. """

    @abc.abstractmethod
    def record_resume_state(self, annotation_run, dump_filename):
        """ Record on the run whatever the import lane needs to resume upload-only off this attempt's
            finished output.

            The reclaim calls this once tool_finished() says the output is complete, because the run's own
            filename fields are only written at the runner's final save - so a run reclaimed before then
            has usable output on disk that its row cannot name. """


class PipelineDef:
    """ Everything the scheduler, dispatcher and upload pipeline need to know about one pipeline type,
        without knowing what its tool is. """

    def __init__(self, runner: AnnotationPipelineRunner,
                 depends_on: Optional[VariantAnnotationPipelineType] = None,
                 blocks_vcf_import: bool = True,
                 enabled_setting: Optional[str] = None):
        self.runner = runner
        # Must be FINISHED on the same AnnotationRangeLock before this may launch. Supplementary
        # pipelines update the VariantAnnotation rows a VEP pipeline wrote, so they queue behind it.
        self.depends_on = depends_on
        # Whether an unfinished run of this type holds up a VCF import - see UploadedVCF.is_fully_annotated
        # (#1656). Supplementary annotation is additive: the VCF is usable without it.
        self.blocks_vcf_import = blocks_vcf_import
        # Settings flag gating whether the scheduler creates runs of this type at all.
        self.enabled_setting = enabled_setting

    @property
    def pipeline_type(self) -> VariantAnnotationPipelineType:
        return self.runner.pipeline_type

    @property
    def enabled(self) -> bool:
        return self.enabled_setting is None or bool(getattr(settings, self.enabled_setting))

    def __str__(self):
        return f"PipelineDef({VariantAnnotationPipelineType(self.pipeline_type).label})"
