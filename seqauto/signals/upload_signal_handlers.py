"""
seqauto's receivers for upload's signals: an imported VCF is linked to the sequencing VCF it was sent for
(seqauto/upload_links.py:link_uploaded_vcf_to_sequencing) and announced again with backend_vcf_import_success_signal
when its import finishes, and a gene coverage file sent for a QCGeneCoverage loads with that record's build and kit.
"""
import logging
import os
from typing import Optional

from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from django.dispatch import receiver

from seqauto.models.models_seqauto import QCGeneCoverage
from seqauto.signals.signals_list import backend_vcf_import_success_signal
from seqauto.upload_links import get_sequencing_record_for_path, link_uploaded_vcf_to_sequencing
from upload.signals.signals import (
    GeneCoverageSource,
    gene_coverage_upload_signal,
    uploaded_vcf_configured_signal,
    vcf_import_success_signal,
)


@receiver(uploaded_vcf_configured_signal)
def uploaded_vcf_configured_handler(sender, uploaded_vcf, upload_step=None, **kwargs):
    link_uploaded_vcf_to_sequencing(uploaded_vcf, upload_step)


@receiver(vcf_import_success_signal)
def vcf_import_success_handler(sender, vcf, **kwargs):
    try:
        backend_vcf = vcf.uploadedvcf.backendvcf
    except ObjectDoesNotExist:
        return
    backend_vcf_import_success_signal.send(sender=os.path.basename(__file__), backend_vcf=backend_vcf)


@receiver(gene_coverage_upload_signal)
def gene_coverage_upload_handler(sender, file_upload, **kwargs) -> Optional[GeneCoverageSource]:
    """ Only the API sets a path - SeqAuto itself loads coverage with GeneCoverageCollection.load_from_file """
    if not (settings.SEQAUTO_ENABLED and file_upload.path):
        return None

    qc_gene_coverage = get_sequencing_record_for_path(file_upload, QCGeneCoverage)
    logging.info("Found seqauto QCGeneCoverage sharing path: %s", qc_gene_coverage.path)

    def link_collection(gene_coverage_collection):
        qc_gene_coverage.gene_coverage_collection = gene_coverage_collection
        qc_gene_coverage.save()

    return GeneCoverageSource(genome_build=qc_gene_coverage.qc.genome_build,
                              enrichment_kit=qc_gene_coverage.qc.sequencing_sample.enrichment_kit,
                              link_collection=link_collection)
