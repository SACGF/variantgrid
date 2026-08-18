from annotation.gene_level_annotation import annotate_gene_level_run
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.pipelines.base import AnnotationPipelineRunner


class GeneLevelRunner(AnnotationPipelineRunner):
    """ Gene-level variants (gene fusions). No VCF to dump and no VEP to run - the annotation is computed
        from the variant's own gene identity, and there are few enough rows to write inline rather than
        handing the run to the import lane, so annotate() takes the run straight to FINISHED. """

    pipeline_type = VariantAnnotationPipelineType.GENE_LEVEL

    def annotate(self, annotation_run, lease_heartbeat=None):
        annotate_gene_level_run(annotation_run)

    def import_results(self, annotation_run):
        """ Never reached - annotate() finishes the run, so it never enters the import lane. """

    def get_output_paths(self, annotation_run, dump_filename) -> list[str]:
        return []  # writes nothing to disk

    def tool_finished(self, annotation_run, dump_filename) -> bool:
        return False  # no expensive output for a reclaim to keep

    def record_resume_state(self, annotation_run, dump_filename):
        """ Never reached - tool_finished() is always False, so a reclaim always scrubs back to CREATED. """
