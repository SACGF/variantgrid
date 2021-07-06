from django.test import TestCase

from snpdb.models import GenomeBuild
from upload.models import ModifiedImportedVariant


class TestModifiedImportedVariant(TestCase):
    def test_format_old_variant(self):
        grch37 = GenomeBuild.get_name_or_alias('GRCh37')
        OLD_VARIANT = "NC_000019.9:536068:G/GTCCTCGTCCTTCCGGGACCCGGGGCGCTGGGAGCCTCACG"
        old_variant_formatted = ModifiedImportedVariant.format_old_variant(OLD_VARIANT, grch37)
        ov = old_variant_formatted.pop()
        chrom = ov.split(":", 1)[0]
        old_chrom = OLD_VARIANT.split(":", 1)[0]
        contig = grch37.chrom_contig_mappings[old_chrom]
        self.assertEqual(chrom, contig.name, "contig converted to chrom name")

    def test_format_old_variant_multi(self):
        grch37 = GenomeBuild.get_name_or_alias('GRCh37')
        OLD_VARIANT_MULTI = "19:536068:G/GTCCTCGTCCTTCCGGGACCCGGGGCGCTGGGAGCCTCACG,NC_000019.9:536046:C/CCCGGGGCGCTGGGAGCCTCACGTCCTCGTCCTTCCGGGAC"
        old_variant_formatted = ModifiedImportedVariant.format_old_variant(OLD_VARIANT_MULTI, grch37)
        self.assertEqual(len(old_variant_formatted), 2)

    def test_format_old_variant_non_decomposed_multi(self):
        """ People may upload a VCF which has already been normalized by VT (but not decomposed) so had the INFO
            already set. As it's normalized already VT won't replace it, thus we'll get a multi but with slashes
            not comma separated """
        grch37 = GenomeBuild.get_name_or_alias('GRCh37')
        OLD_VARIANT_MULTI = "5:132240059:CT/CTT/T"  # Already set
        old_variant_formatted = ModifiedImportedVariant.format_old_variant(OLD_VARIANT_MULTI, grch37)
        self.assertEqual(len(old_variant_formatted), 2)


