from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from annotation.fake_annotation import get_fake_annotation_version
from snpdb.clingen_allele import get_clingen_allele
from snpdb.liftover import (
    _batch_alleles,
    _get_build_liftover_dicts,
    _liftover_using_dest_variant_coordinate,
    _liftover_using_existing_contig,
    _liftover_using_source_variant_coordinate,
    _non_standard_contig_error,
)
from snpdb.models import Allele, AlleleConversionTool, GenomeBuild, VariantCoordinate
from snpdb.tasks.liftover_tasks import _allele_id_batches
from snpdb.tests.utils.mock_clingen_api import MockClinGenAlleleRegistryAPI
from snpdb.tests.utils.vcf_testing_utils import create_mock_allele, slowly_create_test_variant


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


class TestLiftoverBatching(TestCase):
    NUM_ALLELES = 5

    @classmethod
    def setUpTestData(cls):
        cls.alleles = [Allele.objects.create() for _ in range(cls.NUM_ALLELES)]

    @override_settings(LIFTOVER_BATCH_SIZE=2)
    def test_batch_alleles_pages_queryset(self):
        allele_qs = Allele.objects.filter(pk__in=[a.pk for a in self.alleles])
        batches = list(_batch_alleles(allele_qs))
        self.assertEqual([len(b) for b in batches], [2, 2, 1])
        # Paging by pk needs to cover every allele exactly once
        batched_pks = [allele.pk for batch in batches for allele in batch]
        self.assertEqual(sorted(batched_pks), sorted(a.pk for a in self.alleles))

    @override_settings(LIFTOVER_BATCH_SIZE=2)
    def test_batch_alleles_handles_list(self):
        batches = list(_batch_alleles(self.alleles))
        self.assertEqual([len(b) for b in batches], [2, 2, 1])

    @override_settings(LIFTOVER_BATCH_SIZE=2)
    def test_allele_id_batches(self):
        allele_qs = Allele.objects.filter(pk__in=[a.pk for a in self.alleles])
        pks = sorted(a.pk for a in self.alleles)
        expected = [(pks[0], pks[1]), (pks[2], pks[3]), (pks[4], pks[4])]
        self.assertEqual(list(_allele_id_batches(allele_qs)), expected)


class TestLiftoverQueries(TestCase):
    @classmethod
    def setUpTestData(cls):
        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            get_fake_annotation_version(genome_build)

    @staticmethod
    def _create_alleles(num_alleles: int, start_position: int) -> list[Allele]:
        grch37 = GenomeBuild.grch37()
        return [create_mock_allele(slowly_create_test_variant("3", start_position + i, 'A', 'T', grch37), grch37)
                for i in range(num_alleles)]

    @override_settings(LIFTOVER_BCFTOOLS_ENABLED=False)
    def test_queries_do_not_scale_with_alleles(self):
        """ The liftover checks used to query per allele (ClinGen failures, VariantAlleles, ...) """
        grch37 = GenomeBuild.grch37()
        grch38 = GenomeBuild.grch38()

        def _num_queries(alleles) -> int:
            with CaptureQueriesContext(connection) as context:
                _get_build_liftover_dicts(alleles, grch37, [grch38])
            return len(context.captured_queries)

        _num_queries(self._create_alleles(1, 1000))  # Warm up cached contigs/HGVS matcher
        one_allele = _num_queries(self._create_alleles(1, 2000))
        five_alleles = _num_queries(self._create_alleles(5, 3000))
        self.assertEqual(one_allele, five_alleles)
