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
        self.assertEquals(chrom, contig.name, "contig converted to chrom name")

    def test_format_old_variant_multi(self):
        grch37 = GenomeBuild.get_name_or_alias('GRCh37')
        OLD_VARIANT_MULTI = "19:536068:G/GTCCTCGTCCTTCCGGGACCCGGGGCGCTGGGAGCCTCACG,NC_000019.9:536046:C/CCCGGGGCGCTGGGAGCCTCACGTCCTCGTCCTTCCGGGAC"
        old_variant_formatted = ModifiedImportedVariant.format_old_variant(OLD_VARIANT_MULTI, grch37)
        self.assertEquals(len(old_variant_formatted), 2)
