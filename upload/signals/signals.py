"""
upload's hooks for the apps that build on it: vcf_import_success_signal when a VCF's import finishes,
uploaded_vcf_configured_signal once a VCF and its samples exist (seqauto links them to sequencing records) and
gene_coverage_upload_signal for what is known about a gene coverage file before it loads.
"""
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import django.dispatch

from genes.models import GeneCoverageCollection
from snpdb.models import GenomeBuild

vcf_import_success_signal = django.dispatch.Signal()  # args: "vcf"

# Sent when the VCF and its samples are made from the header, and again on a reload, before the extraction upload
# metadata is applied. A receiver raising fails the step. args: "uploaded_vcf", "upload_step"
uploaded_vcf_configured_signal = django.dispatch.Signal()

# Sent before a gene coverage file loads. A receiver that knows the file returns a GeneCoverageSource. args: "file_upload"
gene_coverage_upload_signal = django.dispatch.Signal()


@dataclass
class GeneCoverageSource:
    genome_build: GenomeBuild
    enrichment_kit: Any  # passed to GeneCoverageCollection.load_from_file
    link_collection: Callable[[GeneCoverageCollection], None]
