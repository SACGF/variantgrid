import os
from typing import List

from library.django_utils.django_file_utils import get_import_processing_dir
from library.utils import full_class_name
from library.vcf_utils import write_vcf_from_tuples
from snpdb.models import Variant
from snpdb.models.models_variant import VariantCoordinate
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadedFile, UploadPipeline, UploadStep, \
    UploadedClassificationImport
from upload.models.models_enums import UploadedFileTypes, UploadStepOrigin, \
    UploadStepTaskType, VCFPipelineStage
from upload.upload_processing import process_upload_pipeline
from classification.models.classification import ClassificationImport
from classification.tasks.classification_import_process_variants_task import ClassificationImportProcessVariantsTask


# MAX_VCF_FIELD_LENGTH = 131072
MAX_VCF_FIELD_LENGTH = 1000  # while maximum is much larger than this, it indicated a problem


def _is_safe_for_vcf(variant_coordinate: VariantCoordinate) -> bool:
    if not all([variant_coordinate.chrom, variant_coordinate.pos, variant_coordinate.ref, variant_coordinate.alt]):
        return False

    if len(variant_coordinate.ref) > MAX_VCF_FIELD_LENGTH or len(variant_coordinate.alt) > MAX_VCF_FIELD_LENGTH:
        return False
    return True


def process_classification_import(classification_import: ClassificationImport, import_source):
    """ Classifications are submitted via API with evidence fields (eg HGVS) which we need to
        resolve to a single coordinate - ie link to Variant model.
        If the variant is in the database, link to it, otherwise we need to run it through the VCF import
        pipeline to normalise and/or insert it.

        @see https://github.com/SACGF/variantgrid/wiki/Variant-Classification-Import-and-Liftover

        Batch variant classification submissions are broken up into 1 ClassificationImport per GenomeBuild """

    from classification.models import Classification
    variant_pk_lookup = VariantPKLookup(classification_import.genome_build)
    variant_tuples_by_hash = {}
    classifications_by_hash = {}

    no_variant_qs = classification_import.classification_set.filter(variant__isnull=True)
    classification: Classification
    for classification in no_variant_qs:
        try:
            variant_tuple = classification.get_variant_coordinates_from_evidence(update_flags=True)
            if variant_tuple:
                variant_hash = variant_pk_lookup.add(*variant_tuple)
                if variant_hash:
                    variant_tuples_by_hash[variant_hash] = variant_tuple
                    classifications_by_hash[variant_hash] = classification
            else:
                classification.set_variant(variant=None, message='Could not derive variant coordinates', failed=True)
        except ValueError as ve:
            classification.set_variant(variant=None, message=str(ve), failed=True)

    variant_pk_lookup.batch_check()
    unknown_variant_coordinates = variant_pk_lookup.unknown_variant_coordinates

    for variant_hash, variant_pk in variant_pk_lookup.variant_pk_by_hash.items():
        classification = classifications_by_hash[variant_hash]
        if variant_pk:
            classification.set_variant(Variant.objects.get(pk=variant_pk))
        else:
            variant_tuple = variant_tuples_by_hash[variant_hash]
            if not _is_safe_for_vcf(variant_tuple):
                ref_length = len(variant_tuple.ref) if variant_tuple.ref else 0
                classification.set_variant(None, message=f'Could not process via VCF, ref length = {ref_length}', failed=True)
            else:
                unknown_variant_coordinates.append(variant_tuple)

    _classification_upload_pipeline(classification_import, unknown_variant_coordinates, import_source)


def _classification_upload_pipeline(
        classification_import: ClassificationImport,
        unknown_variant_tuples_list: List[VariantCoordinate],
        import_source):
    """ We always run this even with no variants to insert as we need:
        * create Alleles for variants
        * perform liftover to other builds
        * set c_hgvs cache """
    if unknown_variant_tuples_list:
        working_dir = get_import_processing_dir(classification_import.pk, "classification_import")
        vcf_filename = os.path.join(working_dir, "classification_import.vcf")
        write_vcf_from_tuples(vcf_filename, unknown_variant_tuples_list)
    else:
        vcf_filename = None

    uploaded_file = UploadedFile.objects.create(path=vcf_filename,
                                                import_source=import_source,
                                                name='Variants from API',
                                                user=classification_import.user,
                                                file_type=UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY,
                                                visible=False)

    UploadedClassificationImport.objects.create(uploaded_file=uploaded_file,
                                                classification_import=classification_import)
    upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
    _add_post_data_insertion_upload_steps(upload_pipeline)
    process_upload_pipeline(upload_pipeline)


def _add_post_data_insertion_upload_steps(upload_pipeline: UploadPipeline):
    POST_INSERT_CLASSES = [ClassificationImportProcessVariantsTask]

    sort_order = 10
    for clazz in POST_INSERT_CLASSES:
        class_name = full_class_name(clazz)
        UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                  name=clazz.__name__,
                                  origin=UploadStepOrigin.USER_ADDITION,
                                  sort_order=sort_order,
                                  task_type=UploadStepTaskType.CELERY,
                                  pipeline_stage_dependency=VCFPipelineStage.DATA_INSERTION,
                                  script=class_name)
        sort_order += 1
