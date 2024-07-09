import os
from collections import defaultdict
from time import sleep

from django.contrib.auth.models import User
from django.db.models import QuerySet

from classification.models import ImportedAlleleInfo, ImportedAlleleInfoStatus
from classification.models.classification import ClassificationImport
from classification.models.classification_import_run import ClassificationImportRun
from classification.tasks.classification_import_process_variants_task import ClassificationImportProcessVariantsTask
from library.django_utils.django_file_utils import get_import_processing_dir
from library.genomics.vcf_utils import write_vcf_from_variant_coordinates
from library.utils import full_class_name
from snpdb.models import Variant, ImportSource
from snpdb.models.models_variant import VariantCoordinate
from snpdb.variant_pk_lookup import VariantPKLookup, VariantHash
from upload.models import UploadedFile, UploadPipeline, UploadStep, \
    UploadedClassificationImport
from upload.models.models_enums import UploadedFileTypes, UploadStepOrigin, \
    UploadStepTaskType, VCFPipelineStage
from upload.upload_processing import process_upload_pipeline

# MAX_VCF_FIELD_LENGTH = 131072
MAX_VCF_FIELD_LENGTH = 1000  # while maximum is much larger than this, it indicated a problem


def _is_safe_for_vcf(variant_coordinate: VariantCoordinate) -> bool:
    if not all([variant_coordinate.chrom, variant_coordinate.position, variant_coordinate.ref, variant_coordinate.alt]):
        return False

    if len(variant_coordinate.ref) > MAX_VCF_FIELD_LENGTH or len(variant_coordinate.alt) > MAX_VCF_FIELD_LENGTH:
        return False
    return True


def process_classification_import(classification_import: ClassificationImport, import_source: ImportSource):
    """ Classifications are submitted via API with evidence fields (e.g. HGVS) which we need to
        resolve to a single coordinate - ie link to Variant model.
        If the variant is in the database, link to it, otherwise we need to run it through the VCF import
        pipeline to normalise and/or insert it.

        @see https://github.com/SACGF/variantgrid/wiki/Variant-Classification-Import-and-Liftover
        Batch variant classification submissions are broken up into 1 ClassificationImport per GenomeBuild """

    variant_pk_lookup = VariantPKLookup(classification_import.genome_build)
    variant_coordinates_by_hash: dict[VariantHash, VariantCoordinate] = {}
    allele_info_by_hash: dict[VariantHash, list[ImportedAlleleInfo]] = defaultdict(list)

    # WARNING: previously only matched if the allele_info had no matched_variant
    # (not sure we even have to filter on status of processing, as it shouldn't be part of a classificationImport
    # if it doesn't want to be rematched)
    no_variant_qs = classification_import.importedalleleinfo_set.filter(status=ImportedAlleleInfoStatus.PROCESSING)
    allele_info: ImportedAlleleInfo
    for allele_info in no_variant_qs:
        try:
            if variant_coordinate := allele_info.variant_coordinate_obj:
                variant_hash = variant_pk_lookup.add(variant_coordinate)
                if variant_hash:
                    variant_coordinates_by_hash[variant_hash] = variant_coordinate
                    allele_info_by_hash[variant_hash].append(allele_info)
            else:
                allele_info.set_matching_failed(message='Could not derive variant coordinates')
        except ValueError as ve:
            allele_info.set_matching_failed(message=str(ve))

    variant_pk_lookup.batch_check()
    unknown_variant_coordinates = variant_pk_lookup.unknown_variant_coordinates

    for variant_hash, variant_pk in variant_pk_lookup.variant_pk_by_hash.items():
        for allele_info in allele_info_by_hash[variant_hash]:
            if variant_pk:
                allele_info.set_variant_and_save(matched_variant=Variant.objects.get(pk=variant_pk))
            else:
                variant_coordinate = variant_coordinates_by_hash[variant_hash]
                if not _is_safe_for_vcf(variant_coordinate):
                    ref_length = len(variant_coordinate.ref) if variant_coordinate.ref else 0
                    allele_info.set_matching_failed(message=f'Could not process via VCF, ref length = {ref_length}')
                else:
                    unknown_variant_coordinates.append(variant_coordinate)

    _classification_upload_pipeline(classification_import, unknown_variant_coordinates, import_source)


def _classification_upload_pipeline(
        classification_import: ClassificationImport,
        unknown_variant_coordinates: list[VariantCoordinate],
        import_source: ImportSource):
    """ We always run this even with no variants to insert as we need:
        * create Alleles for variants
        * perform liftover to other builds
        * set c_hgvs cache """
    if unknown_variant_coordinates:
        working_dir = get_import_processing_dir(classification_import.pk, "classification_import")
        vcf_filename = os.path.join(working_dir, "classification_import.vcf")
        write_vcf_from_variant_coordinates(vcf_filename, unknown_variant_coordinates)
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


def variant_matching_dry_run(queryset: QuerySet[ImportedAlleleInfo]):
    for allele_info in queryset.iterator():
        allele_info.dirty_check()


def reattempt_variant_matching(user: User, queryset: QuerySet[ImportedAlleleInfo], clear_existing: bool = False):
    """ @:returns (valid_record_count, invalid_record_count) """
    from classification.models.variant_resolver import VariantResolver
    qs: QuerySet[ImportedAlleleInfo] = queryset.order_by('imported_genome_build_patch_version')
    variant_matcher = VariantResolver(user=user)
    queue_count = 0

    ClassificationImportRun.record_classification_import("admin-variant-rematch")
    try:
        for allele_info in qs:
            if clear_existing:
                allele_info.hard_reset_matching_info()
            queued = variant_matcher.queue_resolve(allele_info)
            if queued:
                queue_count += 1
                if queue_count % 100 == 0:
                    ClassificationImportRun.record_classification_import("admin-variant-rematch", queue_count)
                    queue_count = 0

        return variant_matcher.process_queue()
    finally:
        # got to give time for variant matching to complete, a bit hacky
        sleep(10)
        ClassificationImportRun.record_classification_import("admin-variant-rematch", queue_count, is_complete=True)
