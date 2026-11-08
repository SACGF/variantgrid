"""
VCF processors for Manual Variant Entry, whose VCF ID column holds the ManualVariantEntry pk:
ManualVariantEntryLinkingMixin records a CreatedManualVariant for each record once its variant id is known.
BulkManualVariantEntryLinkingVCFProcessor only inserts variants; BulkNoGenotypeManualVariantEntryLinkingVCFProcessor
also writes the cohort genotypes of the single sample a "create as sample" entry is imported as.
"""
import cyvcf2

from annotation.models.models_manual_variant import CreatedManualVariant
from upload.vcf.bulk_minimal_vcf_processor import BulkMinimalVCFProcessor
from upload.vcf.bulk_no_genotype_vcf_processor import BulkNoGenotypeVCFProcessor


class ManualVariantEntryLinkingMixin:
    """ Each record's ManualVariantEntry pk is kept in step with the processor's variant_hashes, so a batch's
        variant ids line up with the entries that made them """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.manual_variant_entry_ids = []

    def process_entry(self, variant: cyvcf2.Variant):
        # Before super() - processing the entry can complete a batch that includes it
        self.manual_variant_entry_ids.append(int(variant.ID))
        return super().process_entry(variant)

    def batch_variant_ids_known(self, variant_hashes, variant_ids):
        super().batch_variant_ids_known(variant_hashes, variant_ids)
        if len(self.manual_variant_entry_ids) != len(variant_ids):
            raise ValueError(f"{len(self.manual_variant_entry_ids)} manual variant entries for "
                             f"{len(variant_ids)} variants")

        created_manual_variants = [CreatedManualVariant(manual_variant_entry_id=mve_id, variant_id=variant_id)
                                   for mve_id, variant_id in zip(self.manual_variant_entry_ids, variant_ids)]
        if created_manual_variants:
            CreatedManualVariant.objects.bulk_create(created_manual_variants)
        self.manual_variant_entry_ids = []


class BulkManualVariantEntryLinkingVCFProcessor(ManualVariantEntryLinkingMixin, BulkMinimalVCFProcessor):
    """ Inserts the variants and links them to their ManualVariantEntry """


class BulkNoGenotypeManualVariantEntryLinkingVCFProcessor(ManualVariantEntryLinkingMixin, BulkNoGenotypeVCFProcessor):
    """ Imports the entries as a sites-only VCF's single sample, linking variants to their ManualVariantEntry """
