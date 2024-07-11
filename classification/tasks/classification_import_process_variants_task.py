from typing import Optional, Any

from django.conf import settings

from classification.models import ImportedAlleleInfo
from classification.models.classification import ClassificationImport
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import GenomeBuild, ImportSource, Variant, VariantCoordinate, Allele
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import ModifiedImportedVariant, UploadStep
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from variantgrid.celery import app


class ClassificationImportProcessVariantsTask(ImportVCFStepTask):
    """ This is run after the VCF import data insertion stage.
        Variants will be in database, at this stage

        BulkMinimalVCFProcessor will have inserted ModifiedImportedVariant if
        any were normalised during import process """

    def process_items(self, upload_step: UploadStep):
        vc_import = upload_step.uploaded_file.uploadedclassificationimport.classification_import
        genome_build = vc_import.genome_build

        self.link_inserted_variants(genome_build, vc_import, upload_step)

        # Make sure classification variants have ClinGen AlleleIDs
        variants_qs = vc_import.get_variants_qs()
        populate_clingen_alleles_for_variants(genome_build, variants_qs)
        # Schedules liftover (runs as separate tasks)
        liftover_classification_import(vc_import, ImportSource.API)
        # bulk_update_cached_c_hgvs call below won't get liftovers as they're in a separate task.
        # It will be called again in Liftover.complete()
        ImportedAlleleInfo.relink_variants(vc_import)

        return 0  # Unknown how many we have to do so was set to 0

    @staticmethod
    def link_inserted_variants(genome_build: GenomeBuild,
                               classification_import: ClassificationImport,
                               upload_step: UploadStep):
        variant_pk_lookup = VariantPKLookup(genome_build)
        variant_coordinates_by_hash: dict[Any, VariantCoordinate] = {}
        allele_info_by_hash: dict[Any, ImportedAlleleInfo] = {}

        # TODO, should we filter on matched_variant__isnull=True, or on status, or not filter at all so we can rematch
        no_variant_qs = classification_import.importedalleleinfo_set.all()  # .filter(matched_variant__isnull=True)
        for allele_info in no_variant_qs:
            try:
                if variant_coordinate := allele_info.variant_coordinate_obj:
                    variant_hash = variant_pk_lookup.add(variant_coordinate)
                    variant_coordinates_by_hash[variant_hash] = variant_coordinate
                    allele_info_by_hash[variant_hash] = allele_info
                else:
                    pass
            except:
                pass
                # if there's no variant_coordinate, record already knows it's in error, shouldn't have been linked

        # # Create a list of variant tuples for classifications that have no variant set
        # no_variant_qs = classification_import.classification_set.filter(variant__isnull=True)
        # for classification in no_variant_qs:
        #     if variant_coordinate := classification.get_variant_coordinates_from_evidence():
        #         variant_hash = variant_pk_lookup.add(*variant_coordinate)
        #         variant_coordinates_by_hash[variant_hash] = variant_coordinate
        #         classifications_by_hash[variant_hash] = classification
        #     else:
        #         # note this shouldn't happen at this step - to get here get_variant_coordinates_from_evidence
        #         # has to have previously returned a proper value
        #         classification.set_variant_failed_matching(message="Could not derive variant coordinates")

        # Look up variant tuples - normalized ones will be in unknown_variant_coordinates and need to
        # lookup via ModifiedImportedVariant
        variant_pk_lookup.batch_check()

        allele_info_variant_and_message: dict[ImportedAlleleInfo, tuple[Optional[Variant], Optional[str]]] = {}
        for variant_hash, variant_pk in variant_pk_lookup.variant_pk_by_hash.items():
            allele_info = allele_info_by_hash[variant_hash]
            variant = Variant.objects.get(pk=variant_pk)
            allele_info_variant_and_message[allele_info] = (variant, None)

        for variant_coordinate in variant_pk_lookup.unknown_variant_coordinates:
            variant_hash = variant_pk_lookup.get_variant_coordinate_hash(variant_coordinate)
            allele_info = allele_info_by_hash[variant_hash]
            try:
                miv = ModifiedImportedVariant.get_upload_pipeline_unnormalized_variant(upload_step.upload_pipeline,
                                                                                       variant_coordinate)
                variant = miv.variant
                validation_message = f"{miv.old_variant} was normalized to {miv.variant}"
            except ModifiedImportedVariant.DoesNotExist:
                variant_str = " ".join(map(str, variant_coordinate))
                variant = None
                validation_message = f"Variant '{variant_str}' for Allele Info {allele_info.pk} not inserted!"

            allele_info_variant_and_message[allele_info] = (variant, validation_message)

        for allele_info, (variant, message) in allele_info_variant_and_message.items():
            # go via the set method so signals can be called
            if variant:
                allele_info.set_variant_and_save(matched_variant=variant, message=message, force_update=False)
            else:
                allele_info.set_matching_failed(message)


def liftover_classification_import(
        classification_import: ClassificationImport,
        import_source: ImportSource):
    if settings.LIFTOVER_CLASSIFICATIONS:
        variants_qs = classification_import.get_variants_qs()
        allele_qs = Allele.objects.filter(variantallele__variant__in=variants_qs)
        genome_build = classification_import.genome_build
        create_liftover_pipelines(classification_import.user, allele_qs, import_source, genome_build)


ClassificationImportProcessVariantsTask = app.register_task(ClassificationImportProcessVariantsTask())
