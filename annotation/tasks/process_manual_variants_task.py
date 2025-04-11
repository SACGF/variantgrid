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
    """ Variants have already been normalised and inserted via UploadPipeline

        ManualVariantEntry records are made in BulkManualVariantEntryLinkingVCFProcessor
    """

    def process_items(self, upload_step):
        mvec = upload_step.uploaded_file.uploadedmanualvariantentrycollection.collection
        logging.info("ManualVariantsPostInsertTask: mvec_id = %s", mvec)

        if failed_mvec_count := mvec.manualvariantentry_set.filter(createdmanualvariant__isnull=True).count():
            raise ValueError(f"{mvec}: failed to match {failed_mvec_count} to variants")

        mvec.import_status = ImportStatus.SUCCESS
        mvec.save()

        # Ensure ClinGen alleles are set (typically only a small number here,
        # plus people may have created via searching for CAid)
        q = Q(createdmanualvariant__manual_variant_entry__manual_variant_entry_collection=mvec)
        variants_qs = Variant.objects.filter(q)
        populate_clingen_alleles_for_variants(mvec.genome_build, variants_qs)

        return 0  # don't know how many we were supposed to process

ManualVariantsPostInsertTask = app.register_task(ManualVariantsPostInsertTask())
