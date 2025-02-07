from snpdb.clingen_allele import populate_clingen_alleles_for_variants
from snpdb.models import Variant
from upload.vcf.bulk_minimal_vcf_processor import BulkMinimalVCFProcessor


class BulkClinGenAlleleVCFProcessor(BulkMinimalVCFProcessor):
    """ Same as BulkMinimalVCFProcessor but we also insert ClinGen Alleles """

    def batch_handle_variant_ids(self, variant_ids):
        super().batch_handle_variant_ids(variant_ids)

        variants_qs = Variant.objects.filter(pk__in=variant_ids)
        populate_clingen_alleles_for_variants(self.genome_build, variants_qs)
