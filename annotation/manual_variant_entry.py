import os

from django.conf import settings
from django.core.exceptions import PermissionDenied

from annotation.models import ManualVariantEntryType
from annotation.models.models import ManualVariantEntryCollection, ManualVariantEntry
from annotation.tasks.process_manual_variants_task import ManualVariantsPostInsertTask, get_manual_variant_tuples
from library.django_utils.django_file_utils import get_import_processing_dir
from library.utils import full_class_name
from library.vcf_utils import write_vcf_from_tuples
from snpdb.models.models_enums import ImportSource
from snpdb.models.models_genome import GenomeBuild
from upload.models import UploadPipeline, UploadedFile, UploadStep, UploadedManualVariantEntryCollection
from upload.models.models_enums import UploadedFileTypes, UploadStepTaskType, VCFPipelineStage, UploadStepOrigin
from upload.upload_processing import process_upload_pipeline


class CreateManualVariantForbidden(PermissionDenied):
    pass


def can_create_variants(user) -> bool:
    can_create = settings.UPLOAD_ENABLED and settings.VARIANT_MANUAL_CREATE
    can_create &= user.is_superuser or settings.VARIANT_MANUAL_CREATE_BY_NON_ADMIN
    return can_create


def check_can_create_variants(user):
    if not can_create_variants(user):
        raise CreateManualVariantForbidden()


def create_manual_variants(user, genome_build: GenomeBuild, variants_text: str):
    check_can_create_variants(user)
    mvec = ManualVariantEntryCollection.objects.create(user=user,
                                                       genome_build=genome_build)
    variants_list = []
    for i, line in enumerate(variants_text.split('\n')):
        line = line.strip()
        entry_type = ManualVariantEntry.get_entry_type(line)
        kwargs = {"manual_variant_entry_collection": mvec,
                  "line_number": i + 1,
                  "entry_text": line,
                  "entry_type": entry_type}

        if entry_type == ManualVariantEntryType.UNKNOWN:
            kwargs["error_message"] = f"Couldn't determine type of '{line}'"
        mve = ManualVariantEntry.objects.create(**kwargs)

        if entry_type != ManualVariantEntryType.UNKNOWN:
            try:
                variants_list.extend(get_manual_variant_tuples(mve))
            except ValueError as ve:
                mve.error_message = f"Error parsing {entry_type}: '{ve}'"
                mve.save()

    if not variants_list:
        # Pipeline would have just hung forever
        raise ValueError("No valid variants to create!")

    # Because we need to normalise / insert etc, it's easier just to write a VCF
    # and run through upload pipeline
    working_dir = get_import_processing_dir(mvec.pk, "manual_variants")
    vcf_filename = os.path.join(working_dir, "manual_variant_entry.vcf")
    write_vcf_from_tuples(vcf_filename, variants_list)
    uploaded_file = UploadedFile.objects.create(path=vcf_filename,
                                                import_source=ImportSource.WEB,
                                                name='Manual Variant Entry',
                                                user=user,
                                                file_type=UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY,
                                                visible=False)

    UploadedManualVariantEntryCollection.objects.create(uploaded_file=uploaded_file,
                                                        collection=mvec)
    upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)
    add_manual_variant_upload_steps(upload_pipeline)
    process_upload_pipeline(upload_pipeline)
    return mvec


def add_manual_variant_upload_steps(upload_pipeline):

    mv_post_insert_clazz = ManualVariantsPostInsertTask
    class_name = full_class_name(mv_post_insert_clazz)
    UploadStep.objects.create(upload_pipeline=upload_pipeline,
                              name="ManualVariantsPostInsertTask",
                              origin=UploadStepOrigin.USER_ADDITION,
                              sort_order=10,
                              task_type=UploadStepTaskType.CELERY,
                              pipeline_stage_dependency=VCFPipelineStage.DATA_INSERTION,
                              script=class_name)
