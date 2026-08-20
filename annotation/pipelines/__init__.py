"""
The annotation pipelines: what each VariantAnnotationPipelineType actually is.

A pipeline is a tool applied to a class of variant - VEP over short variants, VEP over structural
variants, local computation over gene fusions - so the enum value stored on AnnotationRun is the key
into this registry, and everything that used to be an `if pipeline_type == ...` branch in the tasks,
scheduler and import lives on the PipelineDef or its runner. Adding a tool is adding an entry here.
"""
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.pipelines.annotsv import AnnotSVRunner
from annotation.pipelines.base import AnnotationPipelineRunner, PipelineDef
from annotation.pipelines.gene_level import GeneLevelRunner
from annotation.pipelines.vep import VEPRunner

PIPELINES: dict[VariantAnnotationPipelineType, PipelineDef] = {p.pipeline_type: p for p in [
    PipelineDef(VEPRunner(VariantAnnotationPipelineType.STANDARD)),
    PipelineDef(VEPRunner(VariantAnnotationPipelineType.STRUCTURAL_VARIANT)),
    PipelineDef(GeneLevelRunner(), enabled_setting="ANNOTATION_GENE_LEVEL_ENABLED"),
    PipelineDef(AnnotSVRunner(),
                depends_on=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
                blocks_vcf_import=False,
                enabled_setting="ANNOTATION_ANNOTSV_ENABLED"),
]}


def get_pipeline(pipeline_type) -> PipelineDef:
    return PIPELINES[VariantAnnotationPipelineType(pipeline_type)]


def get_runner(pipeline_type) -> AnnotationPipelineRunner:
    return get_pipeline(pipeline_type).runner


def enabled_pipeline_types() -> list[VariantAnnotationPipelineType]:
    """ Types the scheduler creates runs for on this deployment. """
    return [pt for pt, p in PIPELINES.items() if p.enabled]


def blocking_pipeline_types() -> list[VariantAnnotationPipelineType]:
    """ Types a VCF import waits on - @see UploadedVCF.is_fully_annotated (#1656). A disabled type has
        no runs to wait for, so it can only be a type this deployment actually schedules. """
    return [pt for pt, p in PIPELINES.items() if p.blocks_vcf_import and p.enabled]


def versioned_pipeline_types() -> list[VariantAnnotationPipelineType]:
    """ Types whose runs are scheduled against an AnnotationPipelineVersion (#720) - the ones a tool
        upgrade re-runs, as opposed to VEP's, which roll with the VariantAnnotationVersion. """
    return [pt for pt, p in PIPELINES.items() if p.runner.versioned]


def vep_pipeline_types() -> list[VariantAnnotationPipelineType]:
    """ Types that actually invoke VEP - the only ones get_vep_command means anything for. """
    return [pt for pt, p in PIPELINES.items() if isinstance(p.runner, VEPRunner)]
