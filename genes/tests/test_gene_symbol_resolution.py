"""
Tests for VG's gene-symbol -> transcript resolution path (the code VG owns):
 - DjangoTranscriptDataProvider._get_tags_by_tx_ac() MANE-table fallback (single query)
 - _get_search_hgvs_gene_symbol_transcripts() candidate selection across MANE/all settings,
   which drives VG's provider + override through cdot end-to-end.
"""
from django.test import TestCase, override_settings
from django.test.utils import CaptureQueriesContext
from django.db import connection

from annotation.tests.test_data_fake_genes import _create_fake_gene_version, _insert_transcript_data
from genes.hgvs import HGVSMatcher
from genes.hgvs.biocommons_hgvs.data_provider import DjangoTranscriptDataProvider
from genes.models import MANE, TranscriptVersion
from genes.models_enums import AnnotationConsortium, MANEStatus
from snpdb.models import GenomeBuild
from snpdb.signals.variant_search import _get_search_hgvs_gene_symbol_transcripts


def _make_transcript_version(genome_build, accession, gene_symbol, annotation_consortium,
                             contig, length, tag=None) -> TranscriptVersion:
    gene_id = f"GENE_{accession}"
    gene_version = _create_fake_gene_version(genome_build, gene_id, gene_symbol, annotation_consortium)
    build_data = {
        "url": "fake",
        # single exon [start, end] - get_tx_ac_tags_for_gene ranks by sum(end-start)
        "exons": [[1000, 1000 + length, 0, 1, length, None]],
        "contig": contig,
        "strand": "+",
        "cds_start": 1000,
        "cds_end": 1000 + length,
    }
    if tag is not None:
        build_data["tag"] = tag
    data = {
        "id": accession,
        "cdot": "0.2.27",
        "hgnc": "1101",
        "biotype": ["protein_coding"],
        "gene_name": gene_symbol,
        "start_codon": 1,
        "stop_codon": length,
        "genome_builds": {genome_build.name: build_data},
    }
    return _insert_transcript_data(genome_build, data, gene_version)


class TestGeneSymbolResolution(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.genome_build = GenomeBuild.grch38()
        gb = cls.genome_build
        contig = "17"  # any valid chrom/contig key for this build

        # MANE Select RefSeq transcript - no JSON tag, supplied by the MANE table fallback
        cls.tv_mane_refseq = _make_transcript_version(gb, "NM_000059.4", "BRCA2",
                                                      AnnotationConsortium.REFSEQ, contig, length=3000)
        # Ensembl partner of the MANE Select - also tagged version-insensitively via MANE
        cls.tv_mane_ensembl = _make_transcript_version(gb, "ENST00000380152.8", "BRCA2",
                                                       AnnotationConsortium.ENSEMBL, contig, length=2900)
        # Transcript with a JSON tag (should be preferred over the MANE fallback)
        cls.tv_json_tag = _make_transcript_version(gb, "NM_111111.1", "BRCA2",
                                                   AnnotationConsortium.REFSEQ, contig, length=2000,
                                                   tag="RefSeq_Select,basic")
        # Plain transcript, no tag anywhere
        cls.tv_no_tag = _make_transcript_version(gb, "NM_222222.1", "BRCA2",
                                                 AnnotationConsortium.REFSEQ, contig, length=1000)

        MANE.objects.create(
            symbol=cls.tv_mane_refseq.gene_version.gene_symbol,
            refseq_transcript_version=cls.tv_mane_refseq,
            ensembl_transcript_version=cls.tv_mane_ensembl,
            status=MANEStatus.MANE_SELECT,
        )

    def _bare_data_provider(self) -> DjangoTranscriptDataProvider:
        # Build without __init__ so we don't need a genome FASTA - the tag/ranking
        # paths only touch the DB.
        dp = DjangoTranscriptDataProvider.__new__(DjangoTranscriptDataProvider)
        dp.genome_build = self.genome_build
        return dp

    def test_get_tags_by_tx_ac_mane_fallback(self):
        dp = self._bare_data_provider()
        # cdot passes the canonical versioned accession; MANE matching is version-insensitive
        tx_acs = ["NM_000059.4", "ENST00000380152.8", "NM_111111.1", "NM_222222.1"]
        tags_by_ac = dp._get_tags_by_tx_ac(tx_acs, self.genome_build.name)

        # MANE table fallback supplies the cdot tag spelling for both consortia partners
        self.assertEqual(tags_by_ac["NM_000059.4"], ["MANE_Select"])
        self.assertEqual(tags_by_ac["ENST00000380152.8"], ["MANE_Select"])
        # JSON tag is preferred (not overwritten by the MANE fallback)
        self.assertEqual(tags_by_ac["NM_111111.1"], ["RefSeq_Select", "basic"])
        # No tag anywhere
        self.assertEqual(tags_by_ac["NM_222222.1"], [])

    def test_get_tags_by_tx_ac_single_mane_query(self):
        dp = self._bare_data_provider()
        tx_acs = ["NM_000059.4", "ENST00000380152.8", "NM_111111.1", "NM_222222.1"]
        mane_table = MANE._meta.db_table
        with CaptureQueriesContext(connection) as ctx:
            dp._get_tags_by_tx_ac(tx_acs, self.genome_build.name)
        mane_queries = [q for q in ctx.captured_queries if mane_table in q["sql"]]
        self.assertEqual(len(mane_queries), 1,
                         f"Expected a single MANE query, got {len(mane_queries)}")

    def _matcher_with_bare_dp(self) -> HGVSMatcher:
        # HGVSMatcher whose gene-symbol data provider is the bare (no-FASTA) provider,
        # so rank_gene_symbol_transcripts works without constructing a real converter.
        m = HGVSMatcher.__new__(HGVSMatcher)
        m.genome_build = self.genome_build
        m.hgvs_converter = type("FakeConverter", (), {"hdp": self._bare_data_provider()})()
        return m

    @override_settings(SEARCH_HGVS_GENE_SYMBOL_USE_MANE=True,
                       SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS=False)
    def test_candidate_selection_mane_only(self):
        candidates, has_non_mane = _get_search_hgvs_gene_symbol_transcripts(
            self._matcher_with_bare_dp(), "BRCA2")
        # Only the two MANE Select partners, longest (RefSeq) first
        self.assertEqual(candidates, [
            ("NM_000059.4", MANEStatus.MANE_SELECT.label),
            ("ENST00000380152.8", MANEStatus.MANE_SELECT.label),
        ])
        # NM_111111 (RefSeq_Select) and NM_222222 (untagged) are non-MANE and excluded
        self.assertTrue(has_non_mane)

    @override_settings(SEARCH_HGVS_GENE_SYMBOL_USE_MANE=True,
                       SEARCH_HGVS_GENE_SYMBOL_USE_ALL_TRANSCRIPTS=True)
    def test_candidate_selection_all_transcripts(self):
        candidates, has_non_mane = _get_search_hgvs_gene_symbol_transcripts(
            self._matcher_with_bare_dp(), "BRCA2")
        accessions = [ac for ac, _status in candidates]
        self.assertEqual(set(accessions),
                         {"NM_000059.4", "ENST00000380152.8", "NM_111111.1", "NM_222222.1"})
        self.assertEqual(accessions[0], "NM_000059.4", "Longest MANE Select first")
        status_by_ac = dict(candidates)
        self.assertEqual(status_by_ac["NM_000059.4"], MANEStatus.MANE_SELECT.label)
        self.assertIsNone(status_by_ac["NM_222222.1"])
        # Nothing excluded, so no non-MANE warning
        self.assertFalse(has_non_mane)
