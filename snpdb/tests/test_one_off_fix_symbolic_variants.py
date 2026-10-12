from django.core.management import call_command
from django.test import TestCase, override_settings

from annotation.fake_data import get_fake_annotation_version
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.bcftools_liftover import bcftools_liftover_variant_coordinate
from snpdb.models import AlleleMergeLog, GenomeBuild, Variant, VariantAllele, VariantCoordinate
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant

_POSITION = 128200125
_DEL_LENGTH = 60


class OneOffFixSymbolicVariantsTest(TestCase):
    """ #1358 / #2109 - stored variants brought into line with VARIANT_SYMBOLIC_ALT_SIZE """

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        contig_sequence = cls.grch37.genome_fasta.fasta['3']
        cls.del_ref = contig_sequence[_POSITION - 1:_POSITION + _DEL_LENGTH].upper()

    def _create_legacy_explicit_del(self) -> Variant:
        """ As stored before the threshold came down to 50 """
        with override_settings(VARIANT_SYMBOLIC_ALT_SIZE=1000):
            return slowly_create_test_variant("3", _POSITION, self.del_ref, self.del_ref[0], self.grch37)

    def test_explicit_converted_in_place(self):
        v = self._create_legacy_explicit_del()
        call_command("one_off_fix_symbolic_variants")
        v.refresh_from_db()
        self.assertEqual((VCFSymbolicAllele.DEL, -_DEL_LENGTH, _POSITION, self.del_ref[0]),
                         (v.alt.seq, v.svlen, v.locus.position, v.locus.ref.seq))

    def test_twin_merged_into_canonical(self):
        explicit = self._create_legacy_explicit_del()
        explicit_allele = create_mock_allele(explicit, self.grch37)
        symbolic = slowly_create_test_variant("3", _POSITION, self.del_ref, self.del_ref[0], self.grch37)
        self.assertEqual(VCFSymbolicAllele.DEL, symbolic.alt.seq)
        symbolic_allele = create_mock_allele(symbolic, self.grch37)

        call_command("one_off_fix_symbolic_variants")
        self.assertFalse(Variant.objects.filter(pk=explicit.pk).exists())
        self.assertEqual(symbolic_allele, VariantAllele.objects.get(variant=symbolic).allele)
        self.assertTrue(AlleleMergeLog.objects.filter(old_allele=explicit_allele, new_allele=symbolic_allele,
                                                      success=True).exists())

    def test_bcftools_lifts_short_symbolic_as_explicit(self):
        short_del = VariantCoordinate(chrom="3", position=_POSITION, ref=self.del_ref[0],
                                      alt=VCFSymbolicAllele.DEL, svlen=-_DEL_LENGTH)
        self.assertEqual(VariantCoordinate(chrom="3", position=_POSITION, ref=self.del_ref, alt=self.del_ref[0]),
                         bcftools_liftover_variant_coordinate(short_del, self.grch37))

        long_del = VariantCoordinate(chrom="3", position=_POSITION, ref=self.del_ref[0],
                                     alt=VCFSymbolicAllele.DEL, svlen=-5000)
        self.assertEqual(long_del, bcftools_liftover_variant_coordinate(long_del, self.grch37))
