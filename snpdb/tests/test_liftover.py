from django.db import connection
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext

from annotation.fake_annotation import get_fake_annotation_version
from library.guardian_utils import admin_bot
from snpdb.clingen_allele import get_clingen_allele
from snpdb.liftover import (
    _batch_alleles,
    _get_build_liftover_dicts,
    _liftover_using_dest_variant_coordinate,
    _liftover_using_existing_contig,
    _liftover_using_source_variant_coordinate,
    _non_standard_contig_error,
    _run_liftover_using_same_contig,
    allele_can_attempt_liftover,
)
from snpdb.models import (
    Allele,
    AlleleConversionTool,
    AlleleLiftover,
    AlleleOrigin,
    GenomeBuild,
    LiftoverRun,
    ProcessingStatus,
    VariantAllele,
    VariantCoordinate,
)
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

    def test_retry_conversion_tools_overrides_failed_set(self):
        """ #1273 - a tool that has already failed on an allele is skipped forever unless explicitly retried """
        grch37 = GenomeBuild.grch37()
        grch38 = GenomeBuild.grch38()
        variant_37 = slowly_create_test_variant("3", 128198980, 'A', 'T', grch37)
        VariantAllele.objects.create(variant=variant_37, genome_build=grch37, allele=self.allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE,
                                     allele_linking_tool=AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
        for conversion_tool in [AlleleConversionTool.CLINGEN_ALLELE_REGISTRY,
                                AlleleConversionTool.BCFTOOLS_LIFTOVER]:
            liftover_run = LiftoverRun.objects.create(user=admin_bot(), conversion_tool=conversion_tool,
                                                      genome_build=grch38)
            AlleleLiftover.objects.create(allele=self.allele, liftover=liftover_run,
                                          status=ProcessingStatus.ERROR)

        _existing, needs_pipeline = _get_build_liftover_dicts([self.allele], grch37, [grch38])
        self.assertEqual(dict(needs_pipeline), {}, "Every tool failed, so nothing is attempted")

        retry_tools = {AlleleConversionTool.CLINGEN_ALLELE_REGISTRY}
        _existing, needs_pipeline = _get_build_liftover_dicts([self.allele], grch37, [grch38],
                                                              retry_conversion_tools=retry_tools)
        tools = needs_pipeline[grch38]
        self.assertNotIn(AlleleConversionTool.BCFTOOLS_LIFTOVER, tools, "Tools not retried stay skipped")
        _allele, variant_coordinate, _error = tools[AlleleConversionTool.CLINGEN_ALLELE_REGISTRY][0]
        self.assertEqual(variant_coordinate, self.expected_vc_38)

    def test_allele_can_attempt_liftover_with_retry(self):
        grch38 = GenomeBuild.grch38()
        for conversion_tool in AlleleConversionTool:
            liftover_run = LiftoverRun.objects.create(user=admin_bot(), conversion_tool=conversion_tool,
                                                      genome_build=grch38)
            AlleleLiftover.objects.create(allele=self.allele, liftover=liftover_run,
                                          status=ProcessingStatus.ERROR)

        self.assertFalse(allele_can_attempt_liftover(self.allele, grch38))
        self.assertTrue(allele_can_attempt_liftover(self.allele, grch38,
                                                    retry_conversion_tools=list(AlleleConversionTool)))

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


class TestFailedLiftoverAlleles(TestCase):
    """ #1273 - which alleles a per-tool retry re-attempts """

    @classmethod
    def setUpTestData(cls):
        grch37 = GenomeBuild.grch37()
        grch38 = GenomeBuild.grch38()
        cls.still_missing = create_mock_allele(slowly_create_test_variant("3", 1000, 'A', 'T', grch37), grch37)
        cls.lifted_over = create_mock_allele(slowly_create_test_variant("3", 2000, 'A', 'T', grch37), grch37)
        VariantAllele.objects.create(variant=slowly_create_test_variant("3", 2001, 'A', 'T', grch38),
                                     genome_build=grch38, allele=cls.lifted_over,
                                     origin=AlleleOrigin.LIFTOVER,
                                     allele_linking_tool=AlleleConversionTool.BCFTOOLS_LIFTOVER)

        liftover_run = LiftoverRun.objects.create(user=admin_bot(),
                                                  conversion_tool=AlleleConversionTool.BCFTOOLS_LIFTOVER,
                                                  genome_build=grch38)
        for allele in [cls.still_missing, cls.lifted_over]:
            AlleleLiftover.objects.create(allele=allele, liftover=liftover_run, status=ProcessingStatus.ERROR)

    def test_failed_liftover_for_build(self):
        grch38 = GenomeBuild.grch38()
        qs = Allele.failed_liftover_for_build(grch38, AlleleConversionTool.BCFTOOLS_LIFTOVER)
        self.assertEqual(list(qs), [self.still_missing], "An allele that has since been lifted over is left alone")

        qs = Allele.failed_liftover_for_build(grch38, AlleleConversionTool.CLINGEN_ALLELE_REGISTRY)
        self.assertEqual(list(qs), [], "Only the tool being retried counts")


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


class TestLiftoverSameContig(TestCase):
    """ A build sharing a contig (MT) may already have imported the variant and given it an Allele of its own """

    @classmethod
    def setUpTestData(cls):
        for genome_build in [GenomeBuild.grch37(), GenomeBuild.grch38()]:
            get_fake_annotation_version(genome_build)
        cls.mt_variant = slowly_create_test_variant("MT", 263, 'A', 'G', GenomeBuild.grch37())

    def test_same_contig_liftover_merges_when_dest_variant_already_linked(self):
        grch37 = GenomeBuild.grch37()
        grch38 = GenomeBuild.grch38()
        source_allele = create_mock_allele(self.mt_variant, grch37)
        dest_allele = create_mock_allele(self.mt_variant, grch38)

        liftover = LiftoverRun.objects.create(user=admin_bot(), genome_build=grch38,
                                              conversion_tool=AlleleConversionTool.SAME_CONTIG)
        _run_liftover_using_same_contig(liftover, [(source_allele, self.mt_variant)])

        variant_allele = VariantAllele.objects.get(variant=self.mt_variant, genome_build=grch38)
        surviving_allele_id = min(source_allele.pk, dest_allele.pk)
        self.assertEqual(variant_allele.allele_id, surviving_allele_id)
        self.assertEqual(AlleleLiftover.objects.get(liftover=liftover).status, ProcessingStatus.SKIPPED)

    def test_same_contig_liftover_links_unclaimed_variant(self):
        grch37 = GenomeBuild.grch37()
        grch38 = GenomeBuild.grch38()
        source_allele = create_mock_allele(self.mt_variant, grch37)

        liftover = LiftoverRun.objects.create(user=admin_bot(), genome_build=grch38,
                                              conversion_tool=AlleleConversionTool.SAME_CONTIG)
        _run_liftover_using_same_contig(liftover, [(source_allele, self.mt_variant)])

        variant_allele = VariantAllele.objects.get(variant=self.mt_variant, genome_build=grch38)
        self.assertEqual(variant_allele.allele, source_allele)
        self.assertEqual(AlleleLiftover.objects.get(liftover=liftover).status, ProcessingStatus.SUCCESS)
