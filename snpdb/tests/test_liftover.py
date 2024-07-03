from unittest import TestCase

from snpdb.clingen_allele import get_clingen_allele
from snpdb.liftover import _liftover_using_existing_contig, _liftover_using_dest_variant_coordinate, \
    _liftover_using_source_variant_coordinate
from snpdb.models import GenomeBuild, AlleleConversionTool, VariantCoordinate
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class TestLiftover(TestCase):
    """ Test for library.utils - Django unit tests are run in apps so this needs to be here """

    def setUp(self):
        super().setUp()
        clingen_api = MockClinGenAlleleRegistryAPI()
        clingen_allele = get_clingen_allele("CA10617208", clingen_api=clingen_api)
        self.allele = clingen_allele.allele
        self.expected_vc_37 = VariantCoordinate(chrom='3', position=128198980, ref='A', alt='T')
        self.expected_vc_38 = VariantCoordinate(chrom='3', position=128480137, ref='A', alt='T')

        slowly_create_test_variant("3", 128198980, 'A', 'T', GenomeBuild.grch37())

    def test_liftover_using_existing_variant(self):
        clingen_api = MockClinGenAlleleRegistryAPI()
        # This is on chr3 - it is created for 37 but not 38
        conversion_tool, variant = _liftover_using_existing_contig(self.allele, GenomeBuild.grch38())
        self.assertIsNone(conversion_tool)
        self.assertIsNone(variant)

        # This is MT so should be able to use existing contig
        clingen_allele = get_clingen_allele("CA337095804", clingen_api=clingen_api)
        mt_allele = clingen_allele.allele
        conversion_tool, variant = _liftover_using_existing_contig(mt_allele, GenomeBuild.grch37())
        self.assertEqual(conversion_tool, AlleleConversionTool.SAME_CONTIG)
        self.assertIsNotNone(variant)

    def test_liftover_using_dest_variant_coordinate(self):
        conversion_tool, variant_coordinate_37 = _liftover_using_dest_variant_coordinate(self.allele, GenomeBuild.grch37())
        self.assertEqual(conversion_tool, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
        self.assertEqual(variant_coordinate_37, self.expected_vc_37)

        conversion_tool, variant_coordinate_38 = _liftover_using_dest_variant_coordinate(self.allele, GenomeBuild.grch38())
        self.assertEqual(conversion_tool, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
        self.assertEqual(variant_coordinate_38, self.expected_vc_38)

    def _liftover_using_source_variant_coordinate(self):
        conversion_tool, variant_coordinate_37 = _liftover_using_source_variant_coordinate(self.allele,
                                                                                           source_genome_build=GenomeBuild.grch37(),
                                                                                           dest_genome_build=GenomeBuild.grch38())
        self.assertEqual(conversion_tool, AlleleConversionTool.BCFTOOLS_LIFTOVER)
        self.assertEqual(variant_coordinate_37, self.expected_vc_37)

        conversion_tool, variant_coordinate_38 = _liftover_using_source_variant_coordinate(self.allele,
                                                                                           source_genome_build=GenomeBuild.grch38(),
                                                                                           dest_genome_build=GenomeBuild.grch37())
        self.assertEqual(conversion_tool, AlleleConversionTool.BCFTOOLS_LIFTOVER)
        self.assertEqual(variant_coordinate_38, self.expected_vc_38)
