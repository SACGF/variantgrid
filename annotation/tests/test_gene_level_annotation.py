"""The GENE_LEVEL annotation pipeline - what makes gene lists and comp-het find a fusion."""
from django.test import TestCase

from annotation.annotation_version_querysets import pipeline_type_variant_q
from annotation.fake_annotation import get_fake_annotation_version
from annotation.gene_level_annotation import annotate_gene_level_run
from annotation.models import (
    AnnotationRangeLock,
    AnnotationStatus,
    AnnotationRun,
    VariantAnnotation,
    VariantAnnotationPipelineType,
    VariantGeneOverlap,
    VariantTranscriptAnnotation,
)
from annotation.tests.test_data_fake_genes import _create_fake_gene_version, _insert_transcript_data
from genes.tests.gene_fusion_test_utils import create_gene_fusion
from genes.models import (
    HGNC,
    GeneSymbol,
    HGNCImport,
    ReleaseGeneSymbol,
    ReleaseGeneSymbolGene,
    ReleaseTranscriptVersion,
)
from genes.models_enums import AnnotationConsortium, HGNCStatus
from classification.models.classification_variant_info_models import ImportedAlleleInfo
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion, Variant
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


def _make_gene(genome_build, release, gene_id, gene_symbol, transcript_id, contig, start):
    """ A gene the release knows about, reachable from its symbol """
    gene_version = _create_fake_gene_version(genome_build, gene_id, gene_symbol, AnnotationConsortium.ENSEMBL)
    data = {
        "id": transcript_id,
        "gene_name": gene_symbol,
        "biotype": [],
        "genome_builds": {
            genome_build.name: {
                "url": "fake",
                "exons": [[start, start + 10_000, 0, 1, 10_001, None]],
                "contig": contig,
                "strand": "+",
                "cds_end": start + 10_000,
                "cds_start": start,
            }
        },
    }
    transcript_version = _insert_transcript_data(genome_build, data, gene_version, release)
    ReleaseTranscriptVersion.objects.get_or_create(release=release, transcript_version=transcript_version)

    # What GeneSymbolMatcher builds - symbol -> genes for this release
    release_gene_symbol, _ = ReleaseGeneSymbol.objects.get_or_create(release=release,
                                                                     gene_symbol_id=gene_symbol)
    ReleaseGeneSymbolGene.objects.get_or_create(release_gene_symbol=release_gene_symbol,
                                                gene_id=gene_id)
    return gene_version.gene


class GeneLevelAnnotationTest(TestCase):

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.vav = cls.annotation_version.variant_annotation_version
        release = cls.vav.gene_annotation_release

        hgnc_import = HGNCImport.objects.create()
        cls.hgnc_ids = {}
        for pk, symbol in [(1697, "CD74"), (10261, "ROS1")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
            cls.hgnc_ids[symbol] = pk

        # Deliberately on different contigs - an inter-chromosomal fusion is the case that breaks
        # anchoring on a real chromosome
        cls.cd74_gene = _make_gene(cls.genome_build, release, "ENSG00000000001", "CD74",
                                   "ENST00000000001.1", "5", 149_784_000)
        cls.ros1_gene = _make_gene(cls.genome_build, release, "ENSG00000000002", "ROS1",
                                   "ENST00000000002.1", "6", 117_645_000)

        cls.gene_fusion = create_gene_fusion("CD74", "ROS1")

    def _run_annotation(self) -> AnnotationRun:
        variant = self.gene_fusion.variant
        range_lock = AnnotationRangeLock.objects.create(version=self.vav,
                                                        min_variant=variant, max_variant=variant,
                                                        count=1)
        annotation_run = AnnotationRun.objects.create(
            annotation_range_lock=range_lock,
            pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL)
        annotate_gene_level_run(annotation_run)
        return annotation_run

    def test_pipeline_type_claims_gene_level_only(self):
        variant_qs = Variant.objects.filter(pipeline_type_variant_q(VariantAnnotationPipelineType.GENE_LEVEL))
        self.assertIn(self.gene_fusion.variant, list(variant_qs))

    def test_vep_pipelines_do_not_claim_gene_level(self):
        """ svlen=0 makes a fusion look symbolic, so both VEP types have to subtract it explicitly """
        for pipeline_type in (VariantAnnotationPipelineType.STANDARD,
                              VariantAnnotationPipelineType.STRUCTURAL_VARIANT):
            variant_qs = Variant.objects.filter(pipeline_type_variant_q(pipeline_type))
            self.assertNotIn(self.gene_fusion.variant, list(variant_qs), pipeline_type)

    def test_writes_overlaps_for_both_partners(self):
        """ What makes a gene list on either side find the fusion """
        annotation_run = self._run_annotation()
        gene_ids = set(VariantGeneOverlap.objects.filter(annotation_run=annotation_run)
                       .values_list("gene_id", flat=True))
        self.assertEqual({self.cd74_gene.pk, self.ros1_gene.pk}, gene_ids)

    def test_gene_list_on_the_three_prime_partner_finds_the_fusion(self):
        """ The plan's "a gene list containing ROS1 finds CD74::ROS1" """
        self._run_annotation()
        q = VariantTranscriptAnnotation.get_overlapping_genes_q(self.vav, [self.ros1_gene.pk])
        self.assertIn(self.gene_fusion.variant, list(Variant.objects.filter(q)))

    def test_writes_annotation_with_symbols(self):
        annotation_run = self._run_annotation()
        variant_annotation = VariantAnnotation.objects.get(version=self.vav,
                                                           variant=self.gene_fusion.variant)
        self.assertEqual("CD74", variant_annotation.symbol)
        self.assertEqual("CD74,ROS1", variant_annotation.overlapping_symbols)
        self.assertIsNone(variant_annotation.vep_skipped_reason,
                          "VEP never saw it, so it wasn't skipped - it's blank the way dbNSFP is on an SV")
        self.assertEqual(annotation_run, variant_annotation.annotation_run)
        self.assertTrue(variant_annotation.is_gene_level_annotation)
        self.assertFalse(variant_annotation.has_conservation)

    def test_run_finishes(self):
        """ There is no import lane for gene-level, so the run completes in annotate_gene_level_run """
        annotation_run = self._run_annotation()
        annotation_run.refresh_from_db()
        self.assertEqual(AnnotationStatus.FINISHED, annotation_run.status)
        self.assertEqual(1, annotation_run.dump_count)
        self.assertEqual(1, annotation_run.annotated_count)

    def test_run_with_no_fusions_in_range_finishes(self):
        """ dump_count=0 is what get_status() reads to finish a run with nothing to do - which is
            almost every GENE_LEVEL run, gene-level variants being rare """
        plain_variant = slowly_create_test_variant("3", 128198000, "G", "C", self.genome_build)
        range_lock = AnnotationRangeLock.objects.create(version=self.vav,
                                                        min_variant=plain_variant,
                                                        max_variant=plain_variant, count=1)
        annotation_run = AnnotationRun.objects.create(
            annotation_range_lock=range_lock,
            pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL)

        self.assertEqual(0, annotate_gene_level_run(annotation_run))
        self.assertEqual(AnnotationStatus.FINISHED, annotation_run.status)
        self.assertEqual(0, annotation_run.dump_count)

    def test_rerunnable(self):
        """ Overlaps are rebuilt per annotation version, so running twice is not an error """
        self._run_annotation()
        annotation_run = self._run_annotation()
        self.assertEqual(1, VariantAnnotation.objects.filter(version=self.vav,
                                                             variant=self.gene_fusion.variant).count())
        self.assertEqual(2, VariantGeneOverlap.objects.filter(version=self.vav,
                                                              variant=self.gene_fusion.variant).count(),
                         "one row per partner, not one per run")
        self.assertEqual(annotation_run, VariantAnnotation.objects.get(
            version=self.vav, variant=self.gene_fusion.variant).annotation_run)


class GeneLevelVariantContainmentTest(TestCase):
    """ The rows of the plan's containment table """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(1014, "BCR"), (76, "ABL1")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")
        cls.variant = create_gene_fusion("BCR", "ABL1").variant

    def test_no_clingen_call(self):
        self.assertFalse(self.variant.can_have_clingen_allele)
        self.assertIn("no coordinate", self.variant.clingen_allele_skip_reason())

    def test_no_c_hgvs(self):
        self.assertFalse(self.variant.can_have_c_hgvs)

    def test_as_external_explicit_raises(self):
        genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        with self.assertRaises(ValueError):
            self.variant.coordinate.as_external_explicit(genome_build)

    def test_coordinate_is_gene_level(self):
        self.assertTrue(self.variant.coordinate.is_gene_level)
        self.assertFalse(self.variant.coordinate.can_be_made_explicit)


class GeneFusionClassificationTest(TestCase):
    """ A lab submitting a gene pair as a classification target - @see ImportedAlleleInfo """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(1014, "BCR"), (76, "ABL1")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")

    def _allele_info(self, imported_c_hgvs: str) -> ImportedAlleleInfo:
        return ImportedAlleleInfo.get_or_create(
            imported_c_hgvs=imported_c_hgvs,
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                self.genome_build))

    def test_fusion_string_resolves_to_a_gene_level_coordinate(self):
        """ The coordinate is what the insert pipeline is handed, exactly as for an HGVS that
            resolved - there is one way a variant enters the database """
        allele_info = self._allele_info("BCR::ABL1")
        variant_coordinate = allele_info.variant_coordinate_obj
        self.assertIsNotNone(variant_coordinate)
        self.assertTrue(variant_coordinate.is_gene_level)
        self.assertIn("BCR::ABL1", allele_info.message)

    def test_gene_fusion_is_read_off_the_matched_variant(self):
        allele_info = self._allele_info("BCR::ABL1")
        self.assertIsNone(allele_info.gene_fusion, "nothing is matched until the pipeline inserts it")

        gene_fusion = create_gene_fusion("BCR", "ABL1")
        allele_info.set_variant_and_save(matched_variant=gene_fusion.variant)
        self.assertEqual(gene_fusion, allele_info.gene_fusion)
        self.assertEqual("BCR::ABL1", allele_info.gene_fusion.canonical_str)

    def test_no_c_hgvs_attempted(self):
        allele_info = self._allele_info("BCR::ABL1")
        gene_fusion = create_gene_fusion("BCR", "ABL1")
        allele_info.set_variant_and_save(matched_variant=gene_fusion.variant)
        resolved = allele_info[self.genome_build]
        self.assertIsNone(resolved.c_hgvs, "a fusion sits on no transcript")
        self.assertIsNone(resolved.error, "and that is not an error")
        self.assertEqual("BCR", resolved.gene_symbol_id)

    def test_hgvs_is_left_alone(self):
        """ Checked directly rather than through get_or_create, so the HGVS resolver isn't invoked """
        for imported_c_hgvs in ["NM_000059.3:c.1234A>G", "NC_000007.13:g.140453136A>T"]:
            allele_info = ImportedAlleleInfo(imported_c_hgvs=imported_c_hgvs)
            self.assertFalse(allele_info.resolve_gene_fusion(), imported_c_hgvs)
