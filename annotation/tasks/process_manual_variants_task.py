import logging

from django.db.models.query_utils import Q

from annotation.models import ManualVariantEntryType
from annotation.models.models import ManualVariantEntry
from genes.hgvs import HGVSMatcher
from snpdb.clingen_allele import populate_clingen_alleles_for_variants, get_clingen_alleles_from_external_code
from snpdb.models import Variant, VariantCoordinate
from snpdb.models.models_enums import ImportStatus, ClinGenAlleleExternalRecordType
from upload.models import SimpleVCFImportInfo
from upload.models.models_enums import VCFImportInfoSeverity
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
        hgvs_matcher = HGVSMatcher(mve.genome_build)
        hgvs_string, search_messages = hgvs_matcher.clean_hgvs(mve.entry_text)
        if search_messages:
            mve.warning_message = "\n".join(search_messages)
            mve.save()
        vc = hgvs_matcher.get_variant_coordinate(hgvs_string)
        variant_coordinates.append(vc)
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

        failed_entries = mvec.manualvariantentry_set.filter(createdmanualvariant__isnull=True)
        if failed_mvec_count := failed_entries.count():
            failed_details = []
            for mve in failed_entries:
                detail = f"Line {mve.line_number}: '{mve.entry_text}'"
                if mve.error_message:
                    detail += f" ({mve.error_message})"
                failed_details.append(detail)
            details_str = "; ".join(failed_details)
            message = f"Failed to match {failed_mvec_count} manual variant entries to variants. {details_str}"
            logging.warning("%s: %s", mvec, message)
            SimpleVCFImportInfo.objects.create(
                severity=VCFImportInfoSeverity.WARNING,
                upload_step=upload_step,
                message_string=message,
            )

        mvec.import_status = ImportStatus.SUCCESS
        mvec.save()

        # Ensure ClinGen alleles are set (typically only a small number here,
        # plus people may have created via searching for CAid)
        q = Q(createdmanualvariant__manual_variant_entry__manual_variant_entry_collection=mvec)
        variants_qs = Variant.objects.filter(q)
        populate_clingen_alleles_for_variants(mvec.genome_build, variants_qs)

        return 0  # don't know how many we were supposed to process

ManualVariantsPostInsertTask = app.register_task(ManualVariantsPostInsertTask())
