from django.conf import settings
import cyvcf2

from library.git import Git
from upload.models import VCFImporter
from upload.vcf.abstract_bulk_vcf_processor import AbstractBulkVCFProcessor


class BulkMinimalVCFProcessor(AbstractBulkVCFProcessor):
    """ We have to know the max variant ID in the VCF *even* if we didn't insert
        any variants (as they may have just been inserted and still annotating

        Also inserts Modified Imported Variants
    """

    @staticmethod
    def get_vcf_importer_version():
        """ Get identifier for this importer """

        vcf_importer, _ = VCFImporter.objects.get_or_create(name="MinimalVCFImporter",
                                                            version=1,
                                                            vcf_parser="cyvcf2",
                                                            vcf_parser_version=cyvcf2.__version__,
                                                            code_git_hash=Git(settings.BASE_DIR).hash)
        return vcf_importer

    def process_entry(self, variant):
        ref, alt = self.get_ref_alt(variant)
        variant_hash = self.variant_pk_lookup.get_variant_coordinate_hash(variant.CHROM, variant.POS, ref, alt)
        self.variant_hashes.append(variant_hash)
        self.add_modified_imported_variant(variant, variant_hash)
        self.batch_process_check()
        self.rows_processed += 1

    def batch_handle_variant_ids(self, variant_ids):
        variant_ids_by_hash = dict(zip(self.variant_hashes, variant_ids))
        if self.modified_imported_variants:
            self.process_modified_imported_variants(variant_ids_by_hash)

        self.set_max_variant(self.variant_hashes, variant_ids)
        self.variant_hashes = []

    def batch_process_check(self, minimum_insert_size=None):
        if minimum_insert_size is None:
            minimum_insert_size = self.batch_size

        if self.variant_hashes and len(self.variant_hashes) >= minimum_insert_size:
            variant_ids = self.variant_pk_lookup.get_variant_ids(self.variant_hashes)
            self.batch_handle_variant_ids(variant_ids)

    def finish(self):
        """ This is called at the very end so we can collect any remaining items to process """
        self.batch_process_check(0)  # Insert anything that is there
