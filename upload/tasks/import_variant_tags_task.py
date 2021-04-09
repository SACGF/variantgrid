import logging
import os
import re

import pandas as pd

from analysis.models import VariantTagsImport, ImportedVariantTag, VariantTag, TagLocation
from library.django_utils.django_file_utils import get_import_processing_dir
from library.utils import full_class_name
from library.vcf_utils import write_vcf_from_tuples
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import GenomeBuild, Variant, ImportSource, Tag, VariantAlleleCollectionSource, VariantAllele, \
    VariantAlleleCollectionRecord
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadedVariantTags, UploadedFile, UploadedFileTypes, UploadPipeline, UploadStep, \
    UploadStepOrigin, UploadStepTaskType, VCFPipelineStage
from upload.tasks.import_task import ImportTask
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from upload.upload_processing import process_upload_pipeline
from variantgrid.celery import app


class ImportVariantTagsTask(ImportTask):
    """
        There may be tags with variants we don't know about - so need to insert any unknown

        1. Process the CSV and write ImportedVariantTag (raw string data)
        2. Create a VCF of any variants we don't have.
        3. Insert them as per normal VCF variant only pipeline
        4. Post-process ImportedVariantTag into actual VariantTags
        5. Liftover tags to other builds (VariantAlleleCollectionSource)
    """
    def process_items(self, uploaded_file):
        # There may be variants we don't know about

        df = pd.read_csv(uploaded_file.get_filename())
        NEW_CLASSIFICATION = "New Classification"
        REMOVE_LENGTH = len(NEW_CLASSIFICATION)

        # view_genome_build should all be the same
        view_genome_builds = set(df["view_genome_build"])
        if len(view_genome_builds) != 1:
            raise ValueError(f"Expected 'view_genome_build' column to have 1 unique value, was: {view_genome_builds}")
        genome_build = GenomeBuild.get_name_or_alias(view_genome_builds.pop())
        # Create VariantTagsImport - everything hangs off this
        variant_tags_import = VariantTagsImport.objects.create(user=uploaded_file.user, genome_build=genome_build)
        UploadedVariantTags.objects.create(uploaded_file=uploaded_file, variant_tags_import=variant_tags_import)

        variant_pk_lookup = VariantPKLookup(genome_build)

        imported_tags = []
        for _, row in df.iterrows():
            variant_string = row["variant_string"]
            if variant_string.endswith(NEW_CLASSIFICATION):
                variant_string = variant_string[:-REMOVE_LENGTH].strip()

            variant_tuple = Variant.get_tuple_from_string(variant_string, genome_build=genome_build)
            variant_pk_lookup.add(*variant_tuple)
            node_id = None
            if "node__id" in row:
                node_id = row["node__id"]

            ivt = ImportedVariantTag(variant_tags_import=variant_tags_import,
                                     variant_string=variant_string,
                                     genome_build_string=row["view_genome_build"],
                                     gene_symbol_string=row["variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol"],
                                     tag_string=row["tag__id"],
                                     variant_id=row["variant__id"],
                                     analysis_id=row["analysis__id"],
                                     node_id=node_id,
                                     analysis_name=row["analysis__name"],
                                     user_name=row["user__username"],
                                     created=row["created"])
            imported_tags.append(ivt)

        items_processed = len(imported_tags)
        if imported_tags:
            ImportedVariantTag.objects.bulk_create(imported_tags, batch_size=2000)

        variant_pk_lookup.batch_check()
        if unknown_variant_coordinates := variant_pk_lookup.unknown_variant_coordinates:
            logging.info("There are unknown variants - need to insert, will import tags after")
            self._create_vcf_pipelines(variant_tags_import, unknown_variant_coordinates)
        else:
            logging.info("No unknown variants, can import tags now")
            _create_tags_from_variant_tags_import(variant_tags_import)

        return items_processed

    @staticmethod
    def _create_vcf_pipelines(variant_tags_import, unknown_variant_coordinates):
        # Copy/pasted from classification_import - perhaps refactor?
        working_dir = get_import_processing_dir(variant_tags_import.pk, "variant_tags_import")
        vcf_filename = os.path.join(working_dir, "classification_import.vcf")
        write_vcf_from_tuples(vcf_filename, unknown_variant_coordinates)

        uploaded_file = UploadedFile.objects.create(path=vcf_filename,
                                                    import_source=ImportSource.WEB_UPLOAD,
                                                    name='Variant Tags',
                                                    user=variant_tags_import.user,
                                                    file_type=UploadedFileTypes.VCF_INSERT_VARIANTS_ONLY,
                                                    visible=False)

        upload_pipeline = UploadPipeline.objects.create(uploaded_file=uploaded_file)

        # Once variants are inserted - then actually create the data
        class_name = full_class_name(CreateVariantTagsTask)
        input_filename = f"variant_tags_import={variant_tags_import}"
        UploadStep.objects.create(upload_pipeline=upload_pipeline,
                                  name=class_name,
                                  origin=UploadStepOrigin.USER_ADDITION,
                                  input_filename=input_filename,
                                  sort_order=99,
                                  task_type=UploadStepTaskType.CELERY,
                                  pipeline_stage_dependency=VCFPipelineStage.DATA_INSERTION,
                                  script=class_name)

        process_upload_pipeline(upload_pipeline)


def _create_tags_from_variant_tags_import(variant_tags_import: VariantTagsImport):
    logging.info("_create_tags_from_variant_tags_import: %s!!", variant_tags_import)

    genome_build = variant_tags_import.genome_build

    tag_cache = {}
    variant_tags = []
    for ivt in variant_tags_import.importedvarianttag_set.all():
        tag = tag_cache.get(ivt.tag_string)
        if tag is None:
            tag, _ = Tag.objects.get_or_create(pk=ivt.tag_string)
            tag_cache[tag.pk] = tag

        variant = None
        analysis = None
        node = None
        user = variant_tags_import.user

        # TODO: Attempt to match analysis/node/user
        #ivt.analysis_id
        #ivt.analysis_name
        #ivt.node_id
        #ivt.user_name

        # TODO: We should also look at not creating dupes somehow??

        vt = VariantTag(variant=variant,
                        genome_build=genome_build,
                        tag=tag,
                        analysis=analysis,
                        location=TagLocation.EXTERNAL,
                        imported_variant_tag=ivt,
                        node=node,
                        user=user)
        variant_tags.append(vt)

    if variant_tags:
        logging.info("Creating %d variant tags", len(variant_tags))
        VariantTag.objects.bulk_create(variant_tags, batch_size=2000)

    variant_list = [vt.variant for vt in variant_tags]
    logging.info("Creating liftover pipelines")
    populate_clingen_alleles_for_variants(genome_build, variant_list)
    allele_source = VariantAlleleCollectionSource.objects.create(genome_build=genome_build)
    va_collection_records = []
    for va in VariantAllele.objects.filter(variant__in=variant_list):
        va_collection_records.append(VariantAlleleCollectionRecord(collection=allele_source, variant_allele=va))
    VariantAlleleCollectionRecord.objects.bulk_create(va_collection_records, batch_size=2000)
    logging.info(f"{va_collection_records=}")
    create_liftover_pipelines(variant_tags_import.user, allele_source, ImportSource.WEB, genome_build)


class CreateVariantTagsTask(ImportVCFStepTask):
    """ This is run after the VCF import data insertion stage.
        Variants will be in database, and redis at this stage """

    def process_items(self, upload_step: UploadStep):
        if m := re.match(r"variant_tags_import=(\d+)", upload_step.input_filename):
            variant_tags_import_id = m.group(1)
            variant_tags_import = VariantTagsImport.objects.get(pk=variant_tags_import_id)
            _create_tags_from_variant_tags_import(variant_tags_import)
        else:
            raise ValueError(f"Couldn't extract 'variant_tags_import' PK from '{upload_step.input_filename}'")




ImportVariantTagsTask = app.register_task(ImportVariantTagsTask())  # @UndefinedVariable
CreateVariantTagsTask = app.register_task(CreateVariantTagsTask())  # @UndefinedVariable
