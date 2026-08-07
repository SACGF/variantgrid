"""Unit tests for SVGeneOverlapResolver and BulkVEPVCFAnnotationInserter.add_sv_gene_overlaps.

Issue #1271 — long SVs (>10MB) are skipped by VEP, so we resolve gene overlaps locally
from the VariantAnnotationVersion's gene_annotation_release transcripts.
"""
from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import _create_fake_gene_version, _insert_transcript_data
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import (
    BulkVEPVCFAnnotationInserter,
    SVGeneOverlapResolver,
)
from genes.models import ReleaseTranscriptVersion, TranscriptVersion
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild, VariantCoordinate


def _add_release_transcript_version(release, transcript_version: TranscriptVersion):
    ReleaseTranscriptVersion.objects.get_or_create(release=release,
                                                   transcript_version=transcript_version)


def _make_transcript(genome_build, gene_id, gene_symbol, transcript_id, contig, exons, release):
    gene_version = _create_fake_gene_version(genome_build, gene_id, gene_symbol,
                                             AnnotationConsortium.ENSEMBL)
    data = {
        "id": transcript_id,
        "gene_name": gene_symbol,
        "biotype": [],
        "genome_builds": {
            genome_build.name: {
                "url": "fake",
                "exons": exons,
                "contig": contig,
                "strand": "+",
                "cds_end": exons[-1][1],
                "cds_start": exons[0][0],
            }
        },
    }
    return _insert_transcript_data(genome_build, data, gene_version, release)


class SVGeneOverlapResolverTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh38")
        cls.av = get_fake_annotation_version(cls.genome_build)
        cls.vav = cls.av.variant_annotation_version
        release = cls.vav.gene_annotation_release

        # Two transcripts on the same contig at separate positions.
        # Contig "21" matches the build's chrom_contig_mappings (matches Contig.name).
        cls.tv_a = _make_transcript(cls.genome_build, "ENSG00000000001", "FAKE_A",
                                    "ENST00000000001.1", "21",
                                    exons=[[1_000_000, 1_010_000, 0, 1, 10_001, None]],
                                    release=release)
        cls.tv_b = _make_transcript(cls.genome_build, "ENSG00000000002", "FAKE_B",
                                    "ENST00000000002.1", "21",
                                    exons=[[5_000_000, 5_005_000, 0, 1, 5_001, None]],
                                    release=release)
        _add_release_transcript_version(release, cls.tv_a)
        _add_release_transcript_version(release, cls.tv_b)

        cls.contig_name = cls.tv_a.contig.name

    def test_resolver_overlap_spans_two_transcripts(self):
        resolver = SVGeneOverlapResolver(self.vav)
        # 10 MB deletion spanning both transcripts
        vc = VariantCoordinate(chrom=self.contig_name, position=500_000,
                               ref="N", alt="<DEL>", svlen=-10_000_000)
        symbols, gene_ids = resolver.get_overlaps(vc)
        self.assertEqual(symbols, {"FAKE_A", "FAKE_B"})
        self.assertEqual(gene_ids, {"ENSG00000000001", "ENSG00000000002"})

    def test_resolver_no_overlap(self):
        resolver = SVGeneOverlapResolver(self.vav)
        # SV outside both transcripts
        vc = VariantCoordinate(chrom=self.contig_name, position=20_000_000,
                               ref="N", alt="<DEL>", svlen=-100_000)
        symbols, gene_ids = resolver.get_overlaps(vc)
        self.assertEqual(symbols, set())
        self.assertEqual(gene_ids, set())

    def test_resolver_unknown_contig_returns_empty(self):
        resolver = SVGeneOverlapResolver(self.vav)
        vc = VariantCoordinate(chrom="22", position=1, ref="N", alt="<DEL>", svlen=-1_000)
        symbols, gene_ids = resolver.get_overlaps(vc)
        self.assertEqual(symbols, set())
        self.assertEqual(gene_ids, set())

    def _build_inserter_stub(self):
        """ A stub that supports add_sv_gene_overlaps without standing up the full BulkVEPVCFAnnotationInserter
            (which requires a VEP-annotated VCF and column setup). """

        class _Stub:
            pass

        stub = _Stub()
        stub.sv_gene_overlap_resolver = SVGeneOverlapResolver(self.vav)
        stub.constant_data = {"version_id": self.vav.pk, "annotation_run_id": 999}
        stub.variant_gene_overlap_list = []
        return stub

    def test_add_sv_gene_overlaps_populates_symbols_and_queue(self):
        stub = self._build_inserter_stub()
        vc = VariantCoordinate(chrom=self.contig_name, position=500_000,
                               ref="N", alt="<DEL>", svlen=-10_000_000)
        variant_data = {}
        BulkVEPVCFAnnotationInserter.add_sv_gene_overlaps(stub, variant_id=42,
                                                          variant_coordinate=vc,
                                                          variant_data=variant_data)
        self.assertEqual(variant_data["overlapping_symbols"], "FAKE_A,FAKE_B")
        self.assertEqual(len(stub.variant_gene_overlap_list), 2)
        gene_ids = {row["gene_id"] for row in stub.variant_gene_overlap_list}
        self.assertEqual(gene_ids, {"ENSG00000000001", "ENSG00000000002"})
        for row in stub.variant_gene_overlap_list:
            self.assertEqual(row["variant_id"], 42)
            self.assertEqual(row["version_id"], self.vav.pk)
            self.assertEqual(row["annotation_run_id"], 999)

    def test_add_sv_gene_overlaps_no_overlap_leaves_data_empty(self):
        stub = self._build_inserter_stub()
        vc = VariantCoordinate(chrom=self.contig_name, position=20_000_000,
                               ref="N", alt="<DEL>", svlen=-100_000)
        variant_data = {}
        BulkVEPVCFAnnotationInserter.add_sv_gene_overlaps(stub, variant_id=43,
                                                          variant_coordinate=vc,
                                                          variant_data=variant_data)
        self.assertNotIn("overlapping_symbols", variant_data)
        self.assertEqual(stub.variant_gene_overlap_list, [])
