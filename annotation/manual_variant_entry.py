import os
from typing import Optional

from django.conf import settings
from django.core.exceptions import PermissionDenied

from annotation.models import ManualVariantEntryType
from annotation.models.models_manual_variant import ManualVariantEntry, ManualVariantEntryCollection
from annotation.tasks.process_manual_variants_task import (
    ManualVariantsPostInsertTask,
    get_manual_variant_coordinates,
)
from library.django_utils.django_file_utils import get_import_processing_dir
from library.utils import full_class_name
from snpdb.enums import ImportSource
from snpdb.models.models_genome import GenomeBuild
from snpdb.vcf_utils import write_vcf_from_variant_coordinates
from upload.enums import (
    UploadedFileTypes,
    UploadStepOrigin,
    UploadStepTaskType,
    VCFPipelineStage,
)
from upload.models import (
    FileUpload,
    UploadedManualVariantEntryCollection,
    UploadPipeline,
    UploadStep,
)
from upload.upload_metadata import GENOME_BUILD
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


def create_manual_variants(user, genome_build: GenomeBuild, variants_text: str,
                           sample_name: Optional[str] = None) -> ManualVariantEntryCollection:
    """ sample_name imports the entries as a VCF with one sample of that name, so they can be analysed """
    check_can_create_variants(user)
    mvec = ManualVariantEntryCollection.objects.create(user=user,
                                                       genome_build=genome_build)
    variant_coordinates = []
    vcf_ids = []
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
                for vc in get_manual_variant_coordinates(mve):
                    variant_coordinates.append(vc)
                    vcf_ids.append(mve.id)
            except Exception as ve:
                mve.error_message = f"Error parsing {entry_type}: '{ve}'"
                mve.save()

    if not variant_coordinates:
        # Pipeline would have just hung forever
        raise ValueError("No valid variants to create!")

    # Because we need to normalise / insert etc, it's easier just to write a VCF
    # and run through upload pipeline
    working_dir = get_import_processing_dir(mvec.pk, "manual_variants")
    vcf_filename = os.path.join(working_dir, "manual_variant_entry.vcf")
    write_vcf_from_variant_coordinates(vcf_filename, variant_coordinates, vcf_ids=vcf_ids)
    if sample_name:
        # The VCF header has no contigs to detect the build from
        file_upload_kwargs = {"name": sample_name,
                              "file_type": UploadedFileTypes.MANUAL_VARIANT_ENTRY_SAMPLE,
                              "metadata": {GENOME_BUILD: genome_build.pk}}
    else:
        file_upload_kwargs = {"name": "Manual Variant Entry",
                              "file_type": UploadedFileTypes.MANUAL_VARIANT_ENTRY}
    file_upload = FileUpload.objects.create(path=vcf_filename,
                                            import_source=ImportSource.WEB,
                                            user=user,
                                            **file_upload_kwargs)

    UploadedManualVariantEntryCollection.objects.create(file_upload=file_upload,
                                                        collection=mvec)
    upload_pipeline = UploadPipeline.objects.create(file_upload=file_upload)
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
