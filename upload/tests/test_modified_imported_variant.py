from types import SimpleNamespace

from django.test import TestCase

from snpdb.models import GenomeBuild
from upload.models import ModifiedImportedVariant
from upload.vcf.abstract_bulk_vcf_processor import AbstractBulkVCFProcessor


class TestModifiedImportedVariant(TestCase):
    def test_bcftools_format_old_variant(self):
        grch37 = GenomeBuild.get_name_or_alias('GRCh37')
        OLD_VARIANT = "NC_000019.9|536068|G|GTCCTCGTCCTTCCGGGACCCGGGGCGCTGGGAGCCTCACG"
        old_variant_formatted = ModifiedImportedVariant.bcftools_format_old_variant(OLD_VARIANT,
                                                                                    svlen=None, genome_build=grch37)
        ov = old_variant_formatted.pop()
        chrom = ov.split(":", 1)[0]
        old_chrom = OLD_VARIANT.split("|", 1)[0]
        contig = grch37.chrom_contig_mappings[old_chrom]
        self.assertEqual(chrom, contig.name, "contig converted to chrom name")

    def test_bcftools_format_old_variant_multi(self):
        grch37 = GenomeBuild.get_name_or_alias('GRCh37')
        OLD_VARIANT_MULTI_1 = "19|536068|G|GA,GTCCTCGTCCTTCCGGGACCCGGGGCGCTGGGAGCCTCACG|1"
        old_variant_formatted = ModifiedImportedVariant.bcftools_format_old_variant(OLD_VARIANT_MULTI_1,
                                                                                    svlen=None, genome_build=grch37)
        alt = old_variant_formatted[0].rsplit("/", maxsplit=1)[-1]
        self.assertEqual(alt, "GA")

        OLD_VARIANT_MULTI_2 = "19|536068|G|GA,GTCCTCGTCCTTCCGGGACCCGGGGCGCTGGGAGCCTCACG|2"
        old_variant_formatted = ModifiedImportedVariant.bcftools_format_old_variant(OLD_VARIANT_MULTI_2,
                                                                                    svlen=None, genome_build=grch37)
        alt = old_variant_formatted[0].rsplit("/", maxsplit=1)[-1]
        self.assertEqual(alt, "GTCCTCGTCCTTCCGGGACCCGGGGCGCTGGGAGCCTCACG")

    def test_bcftools_format_zero_alt_index_is_rejected(self):
        """ Alt index is 1-based, so 0 would index the last alt rather than the first. """
        grch37 = GenomeBuild.get_name_or_alias('GRCh37')
        # chrom|pos|ref|alt1,alt2|0  — index 0 is invalid in 1-based numbering
        with self.assertRaises(ValueError):
            ModifiedImportedVariant.bcftools_format_old_variant("1|100|A|C,T|0", svlen=None, genome_build=grch37)


class TestAddModifiedImportedVariant(TestCase):
    """ old_variant is set when bcftools moved or trimmed the record, not when it only split a multi-allelic """

    def _old_variant(self, pos: int, ref: str, alt: str, old_rec: str):
        processor = SimpleNamespace(genome_build=GenomeBuild.get_name_or_alias("GRCh37"))
        record = SimpleNamespace(CHROM="1", POS=pos, REF=ref, ALT=[alt],
                                 INFO={ModifiedImportedVariant.BCFTOOLS_OLD_VARIANT_TAG: old_rec})
        miv_list = []
        AbstractBulkVCFProcessor.add_modified_imported_variant(processor, record, "hash", miv_hash_list=[],
                                                               miv_list=miv_list)
        (_operation, _old_multiallelic, old_variant, _old_variant_formatted, _detail), = miv_list
        return old_variant

    def test_trimmed_in_place(self):
        old_rec = "1|100|GATAT|GATATAT,G|1"
        self.assertEqual(self._old_variant(100, "G", "GAT", old_rec), old_rec)

    def test_moved(self):
        old_rec = "1|104|T|TAT"
        self.assertEqual(self._old_variant(100, "G", "GAT", old_rec), old_rec)

    def test_split_only(self):
        self.assertIsNone(self._old_variant(100, "G", "GAT", "1|100|G|GAT,C|1"))
