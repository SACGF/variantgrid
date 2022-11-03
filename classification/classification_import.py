import os
from collections import defaultdict
from time import sleep
from typing import List, Dict, Tuple

from django.contrib.auth.models import User
from django.db.models import QuerySet

from classification.models.classification import ClassificationImport, Classification
from classification.models.classification_import_run import ClassificationImportRun
from classification.tasks.classification_import_process_variants_task import ClassificationImportProcessVariantsTask
from library.django_utils.django_file_utils import get_import_processing_dir
from library.utils import full_class_name
from library.vcf_utils import write_vcf_from_tuples
from snpdb.models import Variant, ImportSource
from snpdb.models.models_variant import VariantCoordinate
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadedFile, UploadPipeline, UploadStep, \
    UploadedClassificationImport
from upload.models.models_enums import UploadedFileTypes, UploadStepOrigin, \
    UploadStepTaskType, VCFPipelineStage
from upload.upload_processing import process_upload_pipeline

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

    variant_pk_lookup = VariantPKLookup.factory(classification_import.genome_build)
    variant_tuples_by_hash = {}
    classifications_by_hash = defaultdict(list)

    no_variant_qs = classification_import.classification_set.filter(variant__isnull=True)
    classification: Classification
    for classification in no_variant_qs:
        try:
            variant_tuple = classification.get_variant_coordinates_from_evidence(update_flags=True)
            if variant_tuple:
                variant_hash = variant_pk_lookup.add(*variant_tuple)
                if variant_hash:
                    variant_tuples_by_hash[variant_hash] = variant_tuple
                    classifications_by_hash[variant_hash].append(classification)
            else:
                classification.set_variant(variant=None, message='Could not derive variant coordinates', failed=True)
        except ValueError as ve:
            classification.set_variant(variant=None, message=str(ve), failed=True)

    variant_pk_lookup.batch_check()
    unknown_variant_coordinates = variant_pk_lookup.unknown_variant_coordinates

    for variant_hash, variant_pk in variant_pk_lookup.variant_pk_by_hash.items():
        for classification in classifications_by_hash[variant_hash]:
            if variant_pk:
                classification.set_variant(Variant.objects.get(pk=variant_pk))
            else:
                variant_tuple = variant_tuples_by_hash[variant_hash]
                if not _is_safe_for_vcf(variant_tuple):
                    ref_length = len(variant_tuple.ref) if variant_tuple.ref else 0
                    classification.set_variant(None, message=f'Could not process via VCF, ref length = {ref_length}',
                                               failed=True)
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
                                                file_type=UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY)

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


def reattempt_variant_matching(user: User, queryset: QuerySet[Classification]) -> Tuple[int, int]:
    """ @:returns (valid_record_count, invalid_record_count) """

    qs: QuerySet[Classification] = queryset.order_by('evidence__genome_build')
    invalid_record_count = 0
    valid_record_count = 0
    valid_this_loop = 0
    imports_by_genome: Dict[int, ClassificationImport] = {}
    max_size = 100

    def process_outstanding():
        # processes the contents of imports_by_genome
        # and then clears out the appropriate variables
        nonlocal valid_this_loop
        nonlocal valid_record_count
        nonlocal imports_by_genome

        if imports_by_genome:
            for vc_import in imports_by_genome.values():
                process_classification_import(vc_import, ImportSource.API)
            ClassificationImportRun.record_classification_import("admin-variant-rematch", valid_this_loop)

            # reset variables and continue
            valid_record_count += valid_this_loop
            valid_this_loop = 0
            imports_by_genome.clear()

    ClassificationImportRun.record_classification_import("admin-variant-rematch")
    try:

        for vc in qs:
            try:
                vc.revalidate(user=user)
                genome_build = vc.get_genome_build()
                if genome_build.pk not in imports_by_genome:
                    imports_by_genome[genome_build.pk] = ClassificationImport.objects.create(user=user,
                                                                                             genome_build=genome_build)
                vc_import = imports_by_genome[genome_build.pk]
                vc.set_variant(variant=None, message='Admin has re-triggered variant matching')
                vc.classification_import = vc_import
                vc.save()
                valid_this_loop += 1

            except BaseException:
                invalid_record_count += 1

            if valid_this_loop >= max_size:
                process_outstanding()
                sleep(2)  # minor pause to stop the variant matcher from being bombarded

        process_outstanding()

    finally:
        # got to give time for variant matching to complete, a bit hacky
        sleep(10)
        ClassificationImportRun.record_classification_import("admin-variant-rematch", is_complete=True)

    return valid_record_count, invalid_record_count
