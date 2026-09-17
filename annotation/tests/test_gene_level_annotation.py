"""The GENE_LEVEL annotation pipeline - what makes gene lists and comp-het find a gene-level event."""
from django.test import TestCase

from annotation.annotation_version_querysets import (
    get_queryset_for_annotation_version,
    pipeline_type_variant_q,
)
from annotation.fake_annotation import get_fake_annotation_version
from annotation.models.damage_enums import PathogenicityImpact
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
from genes.tests.gene_fusion_test_utils import create_gene_fusion, create_gene_fusion_for_ids
from genes.tests.gene_level_test_utils import create_gene_copy_number_event
from genes.models import (
    HGNC,
    GeneCopyNumberEventKind,
    GeneLevelId,
    GeneFusion,
    GeneSymbol,
    HGNCImport,
    ReleaseGeneSymbol,
    ReleaseGeneSymbolGene,
    ReleaseTranscriptVersion,
    Transcript,
    TranscriptVersion,
)
from genes.models_enums import AnnotationConsortium, HGNCStatus
from library.genomics.vcf_enums import VariantClass
from classification.models.classification_variant_info_models import ImportedAlleleInfo
from snpdb.models import GenomeBuild, GenomeBuildPatchVersion, Variant
from snpdb.variant_queries import get_variant_queryset_for_gene_symbol
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
    return gene_version.gene, transcript_version


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
        cls.cd74_gene, cls.cd74_transcript_version = _make_gene(
            cls.genome_build, release, "ENSG00000000001", "CD74", "ENST00000000001.1", "5", 149_784_000)
        cls.ros1_gene, cls.ros1_transcript_version = _make_gene(
            cls.genome_build, release, "ENSG00000000002", "ROS1", "ENST00000000002.1", "6", 117_645_000)

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

    def test_canonical_c_hgvs_falls_back_to_the_representative_annotation(self):
        """ canonical is VEP's flag, and nothing here sets it - so what the detail page loads comes
            from the representative annotation rather than a transcript row """
        self._run_annotation()
        variant = self.gene_fusion.variant
        self.assertIsNone(variant.get_canonical_transcript_annotation(self.genome_build))
        self.assertEqual("CD74::ROS1", variant.get_canonical_c_hgvs(self.genome_build))

    def test_hgvs_g_is_not_regenerated(self):
        """ The matcher fallback has no coordinate to work from, so it must not be reached """
        self._run_annotation()
        self.assertEqual("CD74::ROS1", VariantAnnotation.get_hgvs_g(self.gene_fusion.variant))

    def test_gene_list_on_the_three_prime_partner_finds_the_fusion(self):
        """ The plan's "a gene list containing ROS1 finds CD74::ROS1" """
        self._run_annotation()
        q = VariantTranscriptAnnotation.get_overlapping_genes_q(self.vav, [self.ros1_gene.pk])
        self.assertIn(self.gene_fusion.variant, list(Variant.objects.filter(q)))

    def test_gene_page_queryset_finds_the_fusion_from_either_partner(self):
        """ The gene page joins VariantGeneOverlap itself rather than going through the transcript
            annotation - @see genes.grids.GeneSymbolVariantsGrid """
        self._run_annotation()
        for symbol in ["CD74", "ROS1"]:
            qs = get_variant_queryset_for_gene_symbol(GeneSymbol.objects.get(pk=symbol),
                                                      self.annotation_version)
            self.assertIn(self.gene_fusion.variant, list(qs), symbol)

    def test_writes_annotation_with_symbols(self):
        annotation_run = self._run_annotation()
        variant_annotation = VariantAnnotation.objects.get(version=self.vav,
                                                           variant=self.gene_fusion.variant)
        self.assertEqual("CD74", variant_annotation.symbol)
        self.assertEqual("CD74,ROS1", variant_annotation.overlapping_symbols)
        self.assertEqual("CD74::ROS1", variant_annotation.hgvs_c)
        self.assertEqual("CD74::ROS1", variant_annotation.hgvs_g,
                         "g.HGVS is the column that is never blank - VICC rather than nothing")
        self.assertEqual("CD74::ROS1", variant_annotation.get_hgvs_c_with_symbol(),
                         "not run through the HGVS parser - it already names both genes")
        self.assertIsNone(variant_annotation.vep_skipped_reason,
                          "VEP never saw it, so it wasn't skipped - it's blank the way dbNSFP is on an SV")
        self.assertEqual(annotation_run, variant_annotation.annotation_run)
        self.assertTrue(variant_annotation.is_gene_level_annotation)
        self.assertFalse(variant_annotation.has_conservation)

    def test_run_before_the_upload_step_creates_the_gene_fusion(self):
        """ The scheduler is kicked as soon as the Variants exist, before the upload pipeline's
            GeneFusion step - a run that gets in first still annotates from the Variant """
        variant = self.gene_fusion.variant
        GeneFusion.objects.filter(variant=variant).delete()
        annotation_run = self._run_annotation()
        gene_ids = set(VariantGeneOverlap.objects.filter(annotation_run=annotation_run)
                       .values_list("gene_id", flat=True))
        self.assertEqual({self.cd74_gene.pk, self.ros1_gene.pk}, gene_ids)
        self.assertTrue(GeneFusion.objects.filter(variant=variant).exists())

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
        """ A second run over the same range finds the variants already annotated, so it writes nothing """
        first_run = self._run_annotation()
        self._run_annotation()
        self.assertEqual(1, VariantAnnotation.objects.filter(version=self.vav,
                                                             variant=self.gene_fusion.variant).count())
        self.assertEqual(2, VariantGeneOverlap.objects.filter(version=self.vav,
                                                              variant=self.gene_fusion.variant).count(),
                         "one row per partner, not one per run")
        self.assertEqual(first_run, VariantAnnotation.objects.get(
            version=self.vav, variant=self.gene_fusion.variant).annotation_run)

    def test_rows_go_into_the_version_partition(self):
        """ Written to the base table they inherit from, they are invisible to every annotation
            queryset - which is how the analysis grid lost its annotation columns """
        self._run_annotation()
        annotation_version = self.vav.get_any_annotation_version()
        for klass in (VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap):
            qs = get_queryset_for_annotation_version(klass, annotation_version)
            self.assertTrue(qs.filter(variant=self.gene_fusion.variant).exists(), klass.__name__)

    def test_representative_gene_and_transcript(self):
        """ No VEP 'pick' to inherit, so the anchor gene and its best transcript stand in """
        self._run_annotation()
        variant_annotation = VariantAnnotation.objects.get(version=self.vav,
                                                           variant=self.gene_fusion.variant)
        self.assertEqual(self.cd74_gene.pk, variant_annotation.gene_id)
        self.assertEqual(self.cd74_transcript_version, variant_annotation.transcript_version)
        self.assertEqual(self.cd74_transcript_version.transcript_id, variant_annotation.transcript_id)

    def test_sets_the_columns_vep_always_fills(self):
        """ consequence/impact/variant_class are non-null on every VEP row, and code downstream
            assumes it - eg gene damage counts does `'missense_variant' in consequence` """
        self._run_annotation()
        variant_annotation = VariantAnnotation.objects.get(version=self.vav,
                                                           variant=self.gene_fusion.variant)
        self.assertEqual("gene_fusion", variant_annotation.consequence)
        self.assertEqual(PathogenicityImpact.HIGH, variant_annotation.impact)
        self.assertEqual(VariantClass.GENE_FUSION, variant_annotation.variant_class)
        self.assertEqual(0, variant_annotation.predictions_num_pathogenic)
        self.assertEqual(0, variant_annotation.predictions_num_benign)

    def test_impact_keeps_fusions_through_a_damage_filter(self):
        """ A null impact drops out of PathogenicityImpact.get_q entirely, so a DamageNode would
            silently lose every fusion """
        self._run_annotation()
        q = PathogenicityImpact.get_q(PathogenicityImpact.MODERATE)
        self.assertIn(self.gene_fusion.variant, list(Variant.objects.filter(q)))

    def test_canonical_says_whether_the_pick_is_tagged(self):
        """ Our transcripts carry no MANE/RefSeq Select tag, so the pick is a fallback, not canonical """
        self._run_annotation()
        variant_annotation = VariantAnnotation.objects.get(version=self.vav,
                                                           variant=self.gene_fusion.variant)
        self.assertFalse(variant_annotation.canonical)

    def test_transcript_version_without_cdot_build_data(self):
        """ Pre-cdot rows have no per-build data - they rank as untagged rather than breaking the run """
        legacy_transcript = Transcript.objects.create(identifier="ENST00000000003",
                                                     annotation_consortium=AnnotationConsortium.ENSEMBL)
        legacy = TranscriptVersion.objects.create(transcript=legacy_transcript, version=1,
                                                 gene_version=self.cd74_transcript_version.gene_version,
                                                 genome_build=self.genome_build,
                                                 contig=self.cd74_transcript_version.contig,
                                                 import_source=self.cd74_transcript_version.import_source,
                                                 data={})
        ReleaseTranscriptVersion.objects.create(release=self.vav.gene_annotation_release,
                                                transcript_version=legacy)
        annotation_run = self._run_annotation()
        transcript_version_ids = set(VariantTranscriptAnnotation.objects.filter(annotation_run=annotation_run)
                                     .values_list("transcript_version_id", flat=True))
        self.assertIn(legacy.pk, transcript_version_ids)

    def test_writes_transcript_annotation_for_both_partners(self):
        """ What the grid export swaps in when an analysis knows its enrichment kit's transcripts """
        annotation_run = self._run_annotation()
        vta_qs = VariantTranscriptAnnotation.objects.filter(annotation_run=annotation_run)
        by_transcript_version = {vta.transcript_version_id: vta for vta in vta_qs}
        self.assertEqual({self.cd74_transcript_version.pk, self.ros1_transcript_version.pk},
                         set(by_transcript_version))
        ros1 = by_transcript_version[self.ros1_transcript_version.pk]
        self.assertEqual("ROS1", ros1.symbol)
        self.assertEqual(self.ros1_gene.pk, ros1.gene_id)
        self.assertEqual("CD74::ROS1", ros1.hgvs_c,
                         "the c.HGVS column keeps its value when a kit's transcript is swapped in")


class GeneLevelAnnotationFromBreakpointGenesTest(TestCase):
    """ A side whose symbol no release carries, reachable only through the genes its breakpoint
        landed in - the ACPP::ETV1 case from the plan """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.vav = cls.annotation_version.variant_annotation_version
        release = cls.vav.gene_annotation_release

        hgnc_import = HGNCImport.objects.create()
        for pk, symbol in [(125, "ACP3"), (3490, "ETV1")]:
            GeneSymbol.objects.get_or_create(symbol=symbol)
            HGNC.objects.create(pk=pk, gene_symbol_id=symbol, hgnc_import=hgnc_import,
                                status=HGNCStatus.APPROVED, approved_name=f"{symbol} approved name")

        cls.acp3_gene, _ = _make_gene(cls.genome_build, release, "ENSG00000000010", "ACP3",
                                      "ENST00000000010.1", "3", 132_000_000)
        cls.etv1_gene, _ = _make_gene(cls.genome_build, release, "ENSG00000000011", "ETV1",
                                      "ENST00000000011.1", "7", 13_930_000)

        # A side minted before the rename was understood: the symbol is one no release matched, and
        # the gene is only known because a breakpoint landed in it
        GeneSymbol.objects.get_or_create(symbol="ACPP")
        cls.anchor = GeneLevelId.objects.create(pk=GeneLevelId.CUSTOM_ID_START,
                                                 symbol_str="ACPP", gene_symbol_id="ACPP")
        cls.anchor.genes.add(cls.acp3_gene)
        cls.partner = GeneLevelId.get_or_create_for_symbol("ETV1", "ETV1", HGNC.objects.get(pk=3490))
        cls.gene_fusion = create_gene_fusion_for_ids(cls.anchor, cls.partner)

    def _run_annotation(self) -> AnnotationRun:
        variant = self.gene_fusion.variant
        range_lock = AnnotationRangeLock.objects.create(version=self.vav, min_variant=variant,
                                                        max_variant=variant, count=1)
        annotation_run = AnnotationRun.objects.create(
            annotation_range_lock=range_lock,
            pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL)
        annotate_gene_level_run(annotation_run)
        return annotation_run

    def test_genes_are_used_where_the_symbol_resolves_to_nothing(self):
        annotation_run = self._run_annotation()
        gene_ids = set(VariantGeneOverlap.objects.filter(annotation_run=annotation_run)
                       .values_list("gene_id", flat=True))
        self.assertEqual({self.acp3_gene.pk, self.etv1_gene.pk}, gene_ids)

    def test_gene_page_queryset_finds_the_fusion_from_the_renamed_side(self):
        """ What a gene list node on ACP3 joins - the plan's failing case """
        self._run_annotation()
        qs = get_variant_queryset_for_gene_symbol(GeneSymbol.objects.get(pk="ACP3"),
                                                  self.annotation_version)
        self.assertIn(self.gene_fusion.variant, list(qs))

    def test_a_side_the_release_cannot_place_writes_no_overlap(self):
        """ Nothing to fall back on - no genes, and a symbol no release matched """
        self.anchor.genes.clear()
        annotation_run = self._run_annotation()
        gene_ids = set(VariantGeneOverlap.objects.filter(annotation_run=annotation_run)
                       .values_list("gene_id", flat=True))
        self.assertEqual({self.etv1_gene.pk}, gene_ids)


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

    def test_c_hgvs_column_is_shown(self):
        """ It carries the VICC nomenclature, so the column has a value rather than being hidden """
        self.assertTrue(self.variant.can_have_c_hgvs)

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
            self.assertFalse(allele_info.resolve_gene_level(), imported_c_hgvs)


class GeneLevelCopyNumberAnnotationTest(TestCase):
    """ A whole-gene copy number event - one gene, and the terms VEP would have used for a
        transcript wholly duplicated or wholly lost """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.vav = cls.annotation_version.variant_annotation_version
        release = cls.vav.gene_annotation_release

        hgnc_import = HGNCImport.objects.create()
        GeneSymbol.objects.get_or_create(symbol="EGFR")
        HGNC.objects.create(pk=3236, gene_symbol_id="EGFR", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="epidermal growth factor receptor")
        cls.egfr_gene, cls.egfr_transcript_version = _make_gene(
            cls.genome_build, release, "ENSG00000000010", "EGFR", "ENST00000000010.1", "7", 55_019_000)

        cls.gain = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN)
        cls.loss = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.LOSS)

    def _run_annotation(self) -> AnnotationRun:
        variants = sorted([self.gain.variant, self.loss.variant], key=lambda v: v.pk)
        range_lock = AnnotationRangeLock.objects.create(version=self.vav,
                                                        min_variant=variants[0], max_variant=variants[-1],
                                                        count=len(variants))
        annotation_run = AnnotationRun.objects.create(
            annotation_range_lock=range_lock,
            pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL)
        annotate_gene_level_run(annotation_run)
        return annotation_run

    def test_writes_the_overlap_for_the_gene(self):
        """ What makes a gene list containing EGFR find 'EGFR amplification' """
        annotation_run = self._run_annotation()
        gene_ids = set(VariantGeneOverlap.objects.filter(annotation_run=annotation_run,
                                                         variant=self.gain.variant)
                       .values_list("gene_id", flat=True))
        self.assertEqual({self.egfr_gene.pk}, gene_ids)

    def test_gene_list_on_the_gene_finds_the_event(self):
        self._run_annotation()
        qs = get_variant_queryset_for_gene_symbol(GeneSymbol.objects.get(pk="EGFR"),
                                                  self.annotation_version)
        self.assertIn(self.gain.variant, list(qs))

    def test_consequence_and_variant_class_per_kind(self):
        """ The grid's kind display and the Effect node both read variant_class """
        self._run_annotation()
        expected = {
            self.gain.variant_id: ("transcript_amplification", VariantClass.COPY_NUMBER_GAIN),
            self.loss.variant_id: ("transcript_ablation", VariantClass.COPY_NUMBER_LOSS),
        }
        for variant_id, (consequence, variant_class) in expected.items():
            variant_annotation = VariantAnnotation.objects.get(version=self.vav, variant_id=variant_id)
            self.assertEqual(consequence, variant_annotation.consequence)
            self.assertEqual(variant_class, variant_annotation.variant_class)
            self.assertEqual(PathogenicityImpact.HIGH, variant_annotation.impact)
            self.assertEqual("EGFR", variant_annotation.symbol)

    def test_hgvs_carries_the_canonical_string(self):
        self._run_annotation()
        self.assertEqual("EGFR amplification", VariantAnnotation.get_hgvs_g(self.gain.variant))
        self.assertEqual("EGFR loss", VariantAnnotation.get_hgvs_g(self.loss.variant))


class GeneCopyNumberClassificationTest(TestCase):
    """ A lab submitting 'EGFR amplification' as a classification target - @see ImportedAlleleInfo """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        hgnc_import = HGNCImport.objects.create()
        GeneSymbol.objects.get_or_create(symbol="EGFR")
        HGNC.objects.create(pk=3236, gene_symbol_id="EGFR", hgnc_import=hgnc_import,
                            status=HGNCStatus.APPROVED, approved_name="epidermal growth factor receptor")

    def _allele_info(self, imported_c_hgvs: str) -> ImportedAlleleInfo:
        return ImportedAlleleInfo.get_or_create(
            imported_c_hgvs=imported_c_hgvs,
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(
                self.genome_build))

    def test_the_words_a_lab_writes_all_reach_one_coordinate(self):
        """ 'amp' and 'gain' are input spellings of what is written out as 'amplification' """
        coordinates = set()
        for imported in ["EGFR amplification", "egfr amp", "EGFR Gain"]:
            allele_info = self._allele_info(imported)
            self.assertTrue(allele_info.variant_coordinate_obj.is_gene_level, imported)
            coordinates.add(allele_info.variant_coordinate)
        self.assertEqual(1, len(coordinates), coordinates)

    def test_loss_and_gain_are_different_coordinates(self):
        self.assertNotEqual(self._allele_info("EGFR amplification").variant_coordinate,
                            self._allele_info("EGFR deletion").variant_coordinate)

    def test_event_is_read_off_the_matched_variant(self):
        allele_info = self._allele_info("EGFR amplification")
        self.assertIsNone(allele_info.gene_copy_number_event,
                          "nothing is matched until the pipeline inserts it")

        event = create_gene_copy_number_event("EGFR", GeneCopyNumberEventKind.GAIN)
        allele_info.set_variant_and_save(matched_variant=event.variant)
        self.assertEqual(event, allele_info.gene_copy_number_event)
        resolved = allele_info[self.genome_build]
        self.assertIsNone(resolved.c_hgvs, "a copy number event sits on no transcript")
        self.assertEqual("EGFR", resolved.gene_symbol_id)

    def test_a_gene_we_do_not_know_mints_nothing(self):
        allele_info = ImportedAlleleInfo(imported_c_hgvs="NOTAGENE amplification")
        self.assertFalse(allele_info.resolve_gene_level())
