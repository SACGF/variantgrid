import logging

from django.db.models.query_utils import Q

from annotation.models import ManualVariantEntryType
from annotation.models.models import CreatedManualVariant, ManualVariantEntry
from genes.hgvs import get_hgvs_variant_coordinate
from snpdb.clingen_allele import populate_clingen_alleles_for_variants, get_clingen_alleles_from_external_code
from snpdb.models import Variant, VariantCoordinate
from snpdb.models.models_enums import ImportStatus, ClinGenAlleleExternalRecordType
from snpdb.variant_pk_lookup import VariantPKLookup
from upload.models import ModifiedImportedVariant
from upload.tasks.vcf.import_vcf_step_task import ImportVCFStepTask
from variantgrid.celery import app


def get_manual_variant_coordinates(mve: ManualVariantEntry) -> list[VariantCoordinate]:
    """ List as dbSNP can return multiple values """
    variant_coordinates = []
    if mve.entry_type == ManualVariantEntryType.DBSNP:
        for clingen_allele in get_clingen_alleles_from_external_code(ClinGenAlleleExternalRecordType.DBSNP_ID,
                                                                     mve.entry_text):
            variant_coordinates.append(clingen_allele.get_variant_coordinate(mve.genome_build))
    elif mve.entry_type == ManualVariantEntryType.HGVS:
        variant_coordinates.append(get_hgvs_variant_coordinate(mve.entry_text, mve.genome_build))
    elif mve.entry_type == ManualVariantEntryType.VARIANT:
        variant_coordinates.append(VariantCoordinate.from_string(mve.entry_text, mve.genome_build))
    else:
        raise ValueError(f"Could not convert entry type of {mve.entry_type}")
    return variant_coordinates


class ManualVariantsPostInsertTask(ImportVCFStepTask):
    """ Variants have already been normalised and inserted via UploadPipeline """

    def process_items(self, upload_step):
        items_processed = 0
        mvec = upload_step.uploaded_file.uploadedmanualvariantentrycollection.collection
        logging.info("ManualVariantsPostInsertTask: mvec_id = %s", mvec)
        variant_pk_lookup = VariantPKLookup(mvec.genome_build)
        variant_coordinate_by_hash = {}
        mvec_id_by_variant_hash = {}

        for mve in mvec.manualvariantentry_set.filter(error_message__isnull=True):
            for variant_coordinate in get_manual_variant_coordinates(mve):
                variant_hash = variant_pk_lookup.add(variant_coordinate)
                variant_coordinate_by_hash[variant_hash] = variant_coordinate
                mvec_id_by_variant_hash[variant_hash] = mve.pk

        logging.info("%d variant hashes", len(variant_coordinate_by_hash))
        created_manual_variants: list[CreatedManualVariant] = []
        variant_pk_lookup.batch_check()
        for variant_hash, variant_pk in variant_pk_lookup.variant_pk_by_hash.items():
            variant_coordinate = variant_coordinate_by_hash[variant_hash]
            mvec_id = mvec_id_by_variant_hash[variant_hash]

            if variant_pk is None:  # must have been normalised
                results = list(ModifiedImportedVariant.get_variants_for_unnormalized_variant(variant_coordinate))
                if not results:
                    raise ValueError("variant hash '%s' wasn't inserted as un-normalised variant!", variant_hash)
                else:
                    for v in results:
                        cmv = CreatedManualVariant(manual_variant_entry_id=mvec_id,
                                                   variant=v)
                        created_manual_variants.append(cmv)

            else:
                cmv = CreatedManualVariant(manual_variant_entry_id=mvec_id,
                                           variant_id=variant_pk)
                created_manual_variants.append(cmv)
            items_processed += 1

        if created_manual_variants:
            CreatedManualVariant.objects.bulk_create(created_manual_variants)

        mvec.import_status = ImportStatus.SUCCESS
        mvec.save()

        # Ensure ClinGen alleles are set (typically only a small number here,
        # plus people may have created via searching for CAid)
        q = Q(createdmanualvariant__manual_variant_entry__manual_variant_entry_collection=mvec)
        variants_qs = Variant.objects.filter(q)
        populate_clingen_alleles_for_variants(mvec.genome_build, variants_qs)

        return 0  # don't know how many we were supposed to process

ManualVariantsPostInsertTask = app.register_task(ManualVariantsPostInsertTask())
