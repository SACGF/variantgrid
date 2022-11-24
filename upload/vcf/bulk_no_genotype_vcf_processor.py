"""
When a VCF has no samples / format etc.

No point complicating BulkGenotypeVCFProcessor with lots of if statements etc
"""

import cyvcf2
from django.conf import settings

from library.git import Git
from patients.models_enums import Zygosity
from upload.models import VCFImporter
from upload.vcf.bulk_genotype_vcf_processor import BulkGenotypeVCFProcessor


class BulkNoGenotypeVCFProcessor(BulkGenotypeVCFProcessor):
    VCF_IMPORTER_VERSION = 1  # Change this if you make a major change to the code.
    # Need to distinguish between no entry and 0, can't use None w/postgres command line inserts
    EMPTY_COHORT_GT_DATA = [
        0,  # ref_count
        0,  # het_count
        0,  # hom_count
        0,  # unk_count
        Zygosity.UNKNOWN_ZYGOSITY,   # samples_zygosity
        None,  # samples_allele_depth
        None,  # samples_allele_frequency
        None,  # samples_read_depth
        None,  # samples_genotype_quality
        None,  # samples_phred_likelihood,
        None,  # samples_filters
    ]

    @staticmethod
    def get_vcf_importer_version():
        """ Get identifier for this importer """

        vcf_importer, _ = VCFImporter.objects.get_or_create(name="PythonNoGenotypeKnownVariantsImporter",
                                                            version=BulkNoGenotypeVCFProcessor.VCF_IMPORTER_VERSION,
                                                            vcf_parser="cyvcf2",
                                                            vcf_parser_version=cyvcf2.__version__,
                                                            code_git_hash=Git(settings.BASE_DIR).hash)
        return vcf_importer

    def process_entry(self, variant):
        # Pre-processed by vcf_filter_unknown_contigs so only recognised contigs present
        ref, alt = self.get_ref_alt(variant)
        alt_hash = self.variant_pk_lookup.get_variant_coordinate_hash(variant.CHROM, variant.POS, ref, alt)

        # Don't need to worry about processing all loci - go straight onto variant lists for insert
        self.variant_hashes.append(alt_hash)
        self.variant_filters.append(self.convert_filters(variant.FILTER))
        self.cohort_genotypes.append(self.EMPTY_COHORT_GT_DATA)
        self.variant_gnomad_af.append(variant.INFO.get("AF"))  # gnomAD

        if self.preprocess_vcf_import_info:
            self.add_modified_imported_variant(variant, alt_hash)
        self.rows_processed += 1
        self.batch_process_check()
