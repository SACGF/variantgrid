from django.test import TestCase

from snpdb.models import GenomeBuild
from upload.models import ModifiedImportedVariant


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

