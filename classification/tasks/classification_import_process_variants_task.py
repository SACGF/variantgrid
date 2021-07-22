from django.conf import settings

from library.log_utils import report_exc_info
from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.liftover import create_liftover_pipelines
from snpdb.models import GenomeBuild, ImportSource, Variant
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import ModifiedImportedVariant, UploadStep
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from classification.models.classification import ClassificationImport,\
    ClassificationImportAlleleSource, Classification
from variantgrid.celery import app


class ClassificationImportProcessVariantsTask(ImportVCFStepTask):
    """ This is run after the VCF import data insertion stage.
        Variants will be in database, and redis at this stage

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
        Classification.bulk_update_cached_c_hgvs(vc_import)

        return 0  # Unknown how many we have to do so was set to 0

    @staticmethod
    def link_inserted_variants(genome_build: GenomeBuild,
                               classification_import: ClassificationImport,
                               upload_step: UploadStep):
        variant_pk_lookup = VariantPKLookup.factory(genome_build)
        variant_tuples_by_hash = {}
        classifications_by_hash = {}

        # Create a list of variant tuples for classifications that have no variant set
        no_variant_qs = classification_import.classification_set.filter(variant__isnull=True)
        for classification in no_variant_qs:
            variant_tuple = classification.get_variant_coordinates_from_evidence()
            if variant_tuple:
                variant_hash = variant_pk_lookup.add(*variant_tuple)
                variant_tuples_by_hash[variant_hash] = variant_tuple
                classifications_by_hash[variant_hash] = classification
            else:
                # note this shouldn't happen at this step - to get here get_variant_coordinates_from_evidence
                # has to have previously returned a proper value
                classification.set_variant(None, message="Could not derive variant coordinates", failed=True)

        # Look up variant tuples - if not exists was normalised during import - lookup ModifiedImportedVariant
        variant_pk_lookup.batch_check()
        for variant_hash, variant_pk in variant_pk_lookup.variant_pk_by_hash.items():
            classification = classifications_by_hash[variant_hash]
            variant_tuple = variant_tuples_by_hash[variant_hash]
            try:
                validation_message = None
                if variant_pk is None:
                    # Not in Redis - could have been normalised during import
                    try:
                        miv = ModifiedImportedVariant.get_upload_pipeline_unnormalized_variant(upload_step.upload_pipeline, *variant_tuple)
                        variant_id = miv.variant.pk
                        validation_message = f"{miv.old_variant} was normalized to {miv.variant}"
                    except ModifiedImportedVariant.DoesNotExist:
                        variant_str = " ".join(map(str, variant_tuple))
                        validation_message = f"Variant '{variant_str}' for classification {classification.pk} was not inserted/in Redis!"

                variant = None
                if variant_pk:
                    variant = Variant.objects.get(pk=variant_pk)

                # go via the set method so signals can be called
                classification.set_variant(variant, message=validation_message, failed=not variant)
            except Exception as e:
                report_exc_info()
                classification.set_variant(None, message=f'Unexpected error during matching {str(e)}', failed=True)


def liftover_classification_import(classification_import: ClassificationImport,
                                           import_source: ImportSource):
    if settings.LIFTOVER_CLASSIFICATIONS:
        allele_source = ClassificationImportAlleleSource.objects.create(classification_import=classification_import)
        genome_build = classification_import.genome_build
        create_liftover_pipelines(classification_import.user, allele_source, import_source, genome_build)


ClassificationImportProcessVariantsTask = app.register_task(ClassificationImportProcessVariantsTask())
