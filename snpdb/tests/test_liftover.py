from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.clingen_allele import get_clingen_allele
from snpdb.liftover import (
    _liftover_using_dest_variant_coordinate,
    _liftover_using_existing_contig,
    _liftover_using_source_variant_coordinate,
    _non_standard_contig_error,
)
from snpdb.models import AlleleConversionTool, GenomeBuild, VariantCoordinate
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class TestLiftover(TestCase):

    @classmethod
    def setUpTestData(cls):
        # Need this for HGVSMatcher
        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            get_fake_annotation_version(genome_build)

        clingen_api = MockClinGenAlleleRegistryAPI()
        clingen_allele = get_clingen_allele("CA10617208", clingen_api=clingen_api)
        cls.allele = clingen_allele.allele
        cls.expected_vc_37 = VariantCoordinate(chrom='3', position=128198980, ref='A', alt='T')
        cls.expected_vc_38 = VariantCoordinate(chrom='3', position=128480137, ref='A', alt='T')

        slowly_create_test_variant("3", 128198980, 'A', 'T', GenomeBuild.grch37())
        # Create the MT for 37 (will be able to re-use for 38)
        slowly_create_test_variant("MT", 263, 'A', 'G', GenomeBuild.grch37())

    def test_liftover_using_existing_variant(self):
        clingen_api = MockClinGenAlleleRegistryAPI()
        # This is on chr3 - it is created for 37 but not 38
        conversion_tool, variant = _liftover_using_existing_contig(self.allele, GenomeBuild.grch38())
        self.assertIsNone(conversion_tool)
        self.assertIsNone(variant)

        # MT variant exists in 37 - shares same contig so should be able to re-use for 38
        clingen_allele = get_clingen_allele("CA337095804", clingen_api=clingen_api)
        mt_allele = clingen_allele.allele
        conversion_tool, variant = _liftover_using_existing_contig(mt_allele, GenomeBuild.grch38())
        self.assertEqual(conversion_tool, AlleleConversionTool.SAME_CONTIG)
        self.assertIsNotNone(variant)

    def test_liftover_using_dest_variant_coordinate(self):
        result = list(_liftover_using_dest_variant_coordinate(self.allele, GenomeBuild.grch37()))[0]
        conversion_tool, variant_coordinate_37, _error_message = result
        self.assertEqual(conversion_tool, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
        self.assertEqual(variant_coordinate_37, self.expected_vc_37)

        result = list(_liftover_using_dest_variant_coordinate(self.allele, GenomeBuild.grch38()))[0]
        conversion_tool, variant_coordinate_38, _error_message = result
        self.assertEqual(conversion_tool, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
        self.assertEqual(variant_coordinate_38, self.expected_vc_38)

    def _liftover_using_source_variant_coordinate(self):
        result = list(_liftover_using_source_variant_coordinate(self.allele,
                                                                source_genome_build=GenomeBuild.grch37(),
                                                                dest_genome_build=GenomeBuild.grch38()))[0]
        conversion_tool, variant_coordinate_37, _error_message = result
        self.assertEqual(conversion_tool, AlleleConversionTool.BCFTOOLS_LIFTOVER)
        self.assertEqual(variant_coordinate_37, self.expected_vc_37)

        result = list(_liftover_using_source_variant_coordinate(self.allele,
                                                                source_genome_build=GenomeBuild.grch38(),
                                                                dest_genome_build=GenomeBuild.grch37()))[0]
        conversion_tool, variant_coordinate_38, _error_message = result
        self.assertEqual(conversion_tool, AlleleConversionTool.BCFTOOLS_LIFTOVER)
        self.assertEqual(variant_coordinate_38, self.expected_vc_38)

    def test_standard_contig_written_to_vcf(self):
        self.assertIsNone(_non_standard_contig_error(GenomeBuild.grch37(), self.expected_vc_37))

    def test_unlocalized_scaffold_rejected_before_vcf(self):
        # issue #1197 - these aren't in the reference fasta, and killed the whole liftover run
        grch37 = GenomeBuild.grch37()
        scaffold = grch37.contigs.get(name='HSCHR1_RANDOM_CTG5')
        scaffold_vc = VariantCoordinate(chrom=scaffold.name, position=1000, ref='A', alt='T')
        error = _non_standard_contig_error(grch37, scaffold_vc)
        self.assertIn(scaffold.get_role_display(), error)
