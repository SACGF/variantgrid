import cyvcf2

from annotation.models import CreatedManualVariant
from library.utils import invert_dict
from snpdb.models import VariantCoordinate
from upload.vcf.bulk_minimal_vcf_processor import BulkMinimalVCFProcessor


class BulkManualVariantEntryLinkingVCFProcessor(BulkMinimalVCFProcessor):
    """ Reads VCF with ManualVariantEntry.pk as ID column - then links to it """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.manual_variant_entry_id_by_variant_hash = {}

    @property
    def genome_build(self):
        return self.upload_step.genome_build

    def process_entry(self, variant: cyvcf2.Variant) -> tuple[VariantCoordinate, str]:
        """ :return variant coordinate, variant_hash """
        variant_coordinate, variant_hash = super().process_entry(variant)
        manual_variant_entry_id = int(variant.ID)
        self.manual_variant_entry_id_by_variant_hash[variant_hash] = manual_variant_entry_id
        return variant_coordinate, variant_hash

    def batch_handle_variant_ids(self, variant_ids):
        variant_ids_by_hash = super().batch_handle_variant_ids(variant_ids)
        variant_hash_by_id = invert_dict(variant_ids_by_hash)

        created_manual_variants: list[CreatedManualVariant] = []
        for variant_id in variant_ids:
            variant_hash = variant_hash_by_id[variant_id]
            mvec_id = self.manual_variant_entry_id_by_variant_hash[variant_hash]
            cmv = CreatedManualVariant(manual_variant_entry_id=mvec_id,
                                       variant_id=variant_id)
            created_manual_variants.append(cmv)

        if created_manual_variants:
            CreatedManualVariant.objects.bulk_create(created_manual_variants)

        self.manual_variant_entry_id_by_variant_hash = {}
