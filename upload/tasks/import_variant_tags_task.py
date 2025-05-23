import logging
from collections import defaultdict

import pandas as pd
from dateutil import parser
from django.utils.timezone import make_aware

from analysis.models import VariantTagsImport, ImportedVariantTag, VariantTag, TagLocation
from library.django_utils import UserMatcher
from library.genomics.vcf_utils import write_vcf_from_variant_coordinates
from library.guardian_utils import assign_permission_to_user_and_groups
from library.pandas_utils import df_nan_to_none
from library.utils import invert_dict
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import GenomeBuild, ImportSource, Tag, VariantAllele, VariantCoordinate, Allele
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import UploadedVariantTags, UploadStep, ModifiedImportedVariant, SimpleVCFImportInfo
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from variantgrid.celery import app


class VariantTagsCreateVCFTask(ImportVCFStepTask):
    """ Write a VCF with variants in VariantTags so they can go through normal insert pipeline """
    # These are the original names (JQGrid js export with raw IDs)
    COL_VARIANT_STRING = "variant_string"
    COL_NODE = "node__id"
    COL_CREATED = "created"
    COL_GENOME_BUILD = "view_genome_build"
    COL_TAG = "tag__id"
    COL_GENE = "variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol"
    COL_VARIANT_ID = "variant__id"
    COL_ANALYSIS_ID = "analysis__id"
    COL_ANALYSIS_NAME = "analysis__name"
    COL_USERNAME = "user__username"

    NEW_COLUMNS = {
        COL_VARIANT_STRING: "Variant",
        COL_NODE: "NodeID",
        COL_CREATED: "Created",
        COL_GENOME_BUILD: "Genome Build",
        COL_TAG: "Tag",
        COL_GENE: "Gene",
        COL_VARIANT_ID: "VariantID",
        COL_ANALYSIS_ID: "AnalysisID",
        COL_ANALYSIS_NAME: "Analysis",
        COL_USERNAME: "Username",
    }

    def process_items(self, upload_step):
        upload_pipeline = upload_step.upload_pipeline
        uploaded_file = upload_pipeline.uploaded_file

        df = pd.read_csv(upload_step.input_filename)
        df = df_nan_to_none(df)
        NEW_CLASSIFICATION = "New Classification"
        REMOVE_LENGTH = len(NEW_CLASSIFICATION)

        if self.COL_GENOME_BUILD not in df.columns:
            df = df.rename(columns=invert_dict(self.NEW_COLUMNS))

        # view_genome_build should all be the same
        view_genome_builds = set(df["view_genome_build"])
        if len(view_genome_builds) != 1:
            raise ValueError(f"Expected 'view_genome_build' column to have 1 unique value, was: {view_genome_builds}")
        genome_build = GenomeBuild.get_name_or_alias(view_genome_builds.pop())

        try:
            # We may be reloading - find previous and delete
            uploaded_variant_tags = UploadedVariantTags.objects.get(uploaded_file=uploaded_file)
            uploaded_variant_tags.variant_tags_import.delete()
            uploaded_variant_tags.delete()
        except UploadedVariantTags.DoesNotExist:
            pass  # Ok to create new ones

        # Create VariantTagsImport - everything hangs off this
        variant_tags_import = VariantTagsImport.objects.create(user=uploaded_file.user, genome_build=genome_build)
        UploadedVariantTags.objects.create(uploaded_file=uploaded_file, variant_tags_import=variant_tags_import)

        variant_coordinates = set()
        imported_tags = []
        num_skipped_records = 0
        num_skipped_with_star = 0
        for _, row in df.iterrows():
            variant_string = row["variant_string"]
            if variant_string.endswith(NEW_CLASSIFICATION):
                variant_string = variant_string[:-REMOVE_LENGTH].strip()

            try:
                variant_coordinate = VariantCoordinate.from_string(variant_string)
            except ValueError:
                num_skipped_records += 1
                if "*" in variant_string:
                    num_skipped_with_star += 1
                logging.warning("Could not convert '%s'", variant_string)
                continue
            variant_coordinates.add(variant_coordinate)
            node_id = None
            if "node__id" in row:
                node_id = row["node__id"]

            gene_symbol = row["variant__variantannotation__transcript_version__gene_version__gene_symbol__symbol"]
            created = parser.parse(row["created"])
            created = make_aware(created)
            ivt = ImportedVariantTag(variant_tags_import=variant_tags_import,
                                     variant_string=variant_string,
                                     genome_build_string=row["view_genome_build"],
                                     gene_symbol_string=gene_symbol,
                                     tag_string=row["tag__id"],
                                     variant_id=row["variant__id"],
                                     analysis_id=row["analysis__id"],
                                     node_id=node_id,
                                     analysis_name=row["analysis__name"],
                                     username=row["user__username"],
                                     created=created)
            imported_tags.append(ivt)

        items_processed = len(imported_tags)
        if imported_tags:
            ImportedVariantTag.objects.bulk_create(imported_tags, batch_size=2000)

        if num_skipped_records:
            message_string = f"Skipped {num_skipped_records} records."
            if num_skipped_with_star:
                message_string += f" ({num_skipped_with_star} containing '*') "
            SimpleVCFImportInfo.objects.create(upload_step=upload_step, message_string=message_string)

        write_vcf_from_variant_coordinates(upload_step.output_filename, variant_coordinates)
        return items_processed


class VariantTagsInsertTask(ImportVCFStepTask):
    """ This is run after the VCF import data insertion stage.
        Variants will be in database at this stage """

    def process_items(self, upload_step: UploadStep):
        uploaded_file = upload_step.upload_pipeline.uploaded_file
        uploaded_variant_tags = uploaded_file.uploadedvarianttags
        variant_tags_import = uploaded_variant_tags.variant_tags_import
        logging.info("_create_tags_from_variant_tags_import: %s!!", variant_tags_import)

        genome_build = variant_tags_import.genome_build
        variant_pk_lookup = VariantPKLookup(genome_build)

        tag_cache = {}
        user_matcher = UserMatcher(default_user=variant_tags_import.user)
        variant_tags = []
        created_date = []
        ivt_variant_hashes = {}
        vts_qs = variant_tags_import.importedvarianttag_set.all()
        for i, ivt in enumerate(vts_qs.order_by("variant_string")):
            variant_coordinate = VariantCoordinate.from_string(ivt.variant_string)
            ivt_variant_hashes[ivt] = variant_pk_lookup.get_variant_coordinate_hash(variant_coordinate)

        # Retrieve them all. Some may be None as they were modified/normalized
        variant_hashes = list(set(ivt_variant_hashes.values()))
        variant_ids_inc_null = variant_pk_lookup.get_variant_ids(variant_hashes, validate_not_null=False)
        variant_ids_by_hash = dict(zip(variant_hashes, variant_ids_inc_null))

        ivt_variants = {}
        for ivt, variant_hash in ivt_variant_hashes.items():
            variant_id = variant_ids_by_hash.get(variant_hash)
            if variant_id is None:
                variant_coordinate = VariantCoordinate.from_string(ivt.variant_string)
                try:
                    variant = ModifiedImportedVariant.get_variant_for_unnormalized_variant(upload_step.upload_pipeline,
                                                                                           variant_coordinate)
                except ModifiedImportedVariant.DoesNotExist as mvi:
                    msg = f"Could not find tag variant '{ivt.variant_string}' as Variant or ModifiedImportedVariant"
                    raise ValueError(msg) from mvi
                variant_id = variant.pk
            ivt_variants[ivt] = variant_id

        logging.info("Loaded variants")
        # The Alleles would have been made in BulkClinGenAlleleVCFProcessor
        variant_ids = set(ivt_variants.values())
        # populate_clingen_alleles_for_variants(genome_build, variants)

        va_qs = VariantAllele.objects.filter(variant__in=variant_ids, genome_build=genome_build)
        allele_id_by_variant_id = dict(va_qs.values_list("variant_id", "allele_id"))
        logging.info("Loaded Alleles")

        for ivt, variant_id in ivt_variants.items():
            tag = tag_cache.get(ivt.tag_string)
            if tag is None:
                tag, _ = Tag.objects.get_or_create(pk=ivt.tag_string)
                tag_cache[tag.pk] = tag

            # We're not going to link analysis/nodes - as probably don't match up across systems
            analysis = None
            node = None
            allele_id = allele_id_by_variant_id.get(variant_id)

            # TODO: We should also look at not creating dupes somehow??
            vt = VariantTag(variant_id=variant_id,
                            allele_id=allele_id,
                            genome_build=genome_build,
                            tag=tag,
                            analysis=analysis,
                            location=TagLocation.EXTERNAL,
                            imported_variant_tag=ivt,
                            node=node,
                            user=user_matcher.get_user(ivt.username))

            created_date.append(ivt.created)
            variant_tags.append(vt)

        variant_list = []
        if variant_tags:
            logging.info("Creating %d variant tags", len(variant_tags))
            # VariantTag.created will be set by auto_now_add (no way to stop this)
            variant_tags = VariantTag.objects.bulk_create(variant_tags, batch_size=2000)

            # As we bulk created - permissions weren't done - need to make them
            tags_by_user = defaultdict(list)
            for vt in variant_tags:
                tags_by_user[vt.user].append(vt)

            for user, user_tags in tags_by_user.items():
                assign_permission_to_user_and_groups(user, user_tags)

            # Update date to be from imported date
            for vt, created in zip(variant_tags, created_date):
                vt.created = created
                variant_list.append(vt.variant)
            VariantTag.objects.bulk_update(variant_tags, fields=["created"], batch_size=2000)

        logging.info("Creating liftover pipelines")
        allele_qs = Allele.objects.filter(variantallele__variant__in=variant_list)
        create_liftover_pipelines(variant_tags_import.user, allele_qs, ImportSource.WEB, genome_build)


VariantTagsCreateVCFTask = app.register_task(VariantTagsCreateVCFTask())
VariantTagsInsertTask = app.register_task(VariantTagsInsertTask())
