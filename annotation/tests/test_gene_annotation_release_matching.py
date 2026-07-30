from django.test import TestCase

from annotation.gene_annotation_release_matching import VepGeneSetVersions, match_gene_annotation_release
from genes.models_enums import AnnotationConsortium

NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606"
ENSEMBL = "https://ftp.ensembl.org/pub"

# The real VEP version strings each of these gene sets reports
REFSEQ_GRCH38_RS_2025_08 = "GCF_000001405.40-RS_2025_08 - GCF_000001405.40_GRCh38.p14_genomic.gff"
REFSEQ_GRCH38_RS_2023_10 = "GCF_000001405.40-RS_2023_10 - GCF_000001405.40_GRCh38.p14_genomic.gff"
REFSEQ_GRCH38_110 = "110 - GCF_000001405.40_GRCh38.p14_genomic.gff"
REFSEQ_GRCH37 = "105.20220307 - GCF_000001405.25_GRCh37.p13_genomic.gff"
REFSEQ_GRCH37_LEGACY = "2020-10-26 17:03:42 - GCF_000001405.25_GRCh37.p13_genomic.gff"


def _refseq(genome_build_name, refseq, vep=116):
    return VepGeneSetVersions(genome_build_name=genome_build_name,
                              annotation_consortium=AnnotationConsortium.REFSEQ,
                              refseq=refseq, vep=vep, vep_cache=vep, genebuild="")


def _ensembl(genome_build_name, vep=116, vep_cache=None, genebuild=""):
    return VepGeneSetVersions(genome_build_name=genome_build_name,
                              annotation_consortium=AnnotationConsortium.ENSEMBL,
                              refseq="", vep=vep, vep_cache=vep_cache, genebuild=genebuild)


class ReleaseTokenTest(TestCase):
    def test_refseq_assembly_prefixed_release(self):
        self.assertEqual("RS_2025_08", _refseq("GRCh38", REFSEQ_GRCH38_RS_2025_08).release_token)

    def test_refseq_numbered_release(self):
        self.assertEqual("110", _refseq("GRCh38", REFSEQ_GRCH38_110).release_token)

    def test_refseq_grch37(self):
        self.assertEqual("105.20220307", _refseq("GRCh37", REFSEQ_GRCH37).release_token)

    def test_refseq_grch37_legacy_timestamp(self):
        """ One VEP reports a timestamp rather than an annotation release """
        self.assertEqual("105.20201022", _refseq("GRCh37", REFSEQ_GRCH37_LEGACY).release_token)

    def test_refseq_missing(self):
        self.assertIsNone(_refseq("GRCh38", "").release_token)

    def test_ensembl_grch38_uses_cache_version(self):
        self.assertEqual("115", _ensembl("GRCh38", vep=116, vep_cache=115).release_token)

    def test_ensembl_grch38_falls_back_to_code_version(self):
        self.assertEqual("116", _ensembl("GRCh38", vep=116, vep_cache=None).release_token)

    def test_ensembl_grch37_is_pinned(self):
        self.assertEqual("87", _ensembl("GRCh37", vep=116, vep_cache=116).release_token)

    def test_ensembl_t2t_uses_genebuild(self):
        self.assertEqual("2022_06", _ensembl("T2T-CHM13v2.0", vep=112, vep_cache=107,
                                             genebuild="2022-06").release_token)


class GffUrlTest(TestCase):
    """ These must equal the urls cdot records on its transcripts """

    def test_refseq_assembly_prefixed_release_has_no_assembly_dir(self):
        self.assertEqual(
            f"{NCBI}/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
            _refseq("GRCh38", REFSEQ_GRCH38_RS_2025_08).gff_url)

    def test_refseq_numbered_release_nests_under_assembly_dir(self):
        self.assertEqual(
            f"{NCBI}/110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
            _refseq("GRCh38", REFSEQ_GRCH38_110).gff_url)

    def test_refseq_grch37(self):
        self.assertEqual(
            f"{NCBI}/105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
            _refseq("GRCh37", REFSEQ_GRCH37).gff_url)

    def test_ensembl_grch38_is_gtf(self):
        self.assertEqual(f"{ENSEMBL}/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz",
                         _ensembl("GRCh38", vep=116, vep_cache=116).gff_url)

    def test_ensembl_grch37_is_gtf(self):
        self.assertEqual(f"{ENSEMBL}/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz",
                         _ensembl("GRCh37").gff_url)

    def test_t2t_url_not_derived(self):
        """ Rapid-release urls embed a GCA accession we don't hold - token matching covers these """
        self.assertIsNone(_ensembl("T2T-CHM13v2.0", genebuild="2022-06").gff_url)


class MatchGeneAnnotationReleaseTest(TestCase):
    def test_exact_url_match(self):
        candidates = {
            7: f"{NCBI}/GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
            8: f"{NCBI}/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
        }
        match = match_gene_annotation_release(candidates, _refseq("GRCh38", REFSEQ_GRCH38_RS_2025_08))
        self.assertEqual(8, match.unique_release_id)

    def test_matches_across_url_scheme(self):
        candidates = {8: f"{NCBI}/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
                         .replace("https://", "ftp://")}
        match = match_gene_annotation_release(candidates, _refseq("GRCh38", REFSEQ_GRCH38_RS_2025_08))
        self.assertEqual(8, match.unique_release_id)

    def test_matches_ensembl_gff3_recorded_as_gtf(self):
        candidates = {9: f"{ENSEMBL}/release-116/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz"}
        match = match_gene_annotation_release(candidates, _ensembl("GRCh38", vep=116, vep_cache=116))
        self.assertEqual(9, match.unique_release_id)

    def test_token_fallback_for_t2t(self):
        candidates = {5: ("https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/"
                          "ensembl/geneset/2022_06/Homo_sapiens-GCA_009914755.4-2022_06-genes.gtf.gz"),
                      6: ("https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/"
                          "ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz")}
        match = match_gene_annotation_release(candidates, _ensembl("T2T-CHM13v2.0", vep=112, vep_cache=107,
                                                                  genebuild="2022-06"))
        self.assertEqual(5, match.unique_release_id)

    def test_no_match_when_release_absent(self):
        candidates = {7: f"{NCBI}/GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"}
        match = match_gene_annotation_release(candidates, _refseq("GRCh38", REFSEQ_GRCH38_RS_2025_08))
        self.assertEqual([], match.release_ids)
        self.assertIsNone(match.unique_release_id)

    def test_release_token_does_not_match_a_longer_one(self):
        """ '110' must not match a release directory that merely contains those digits """
        candidates = {7: f"{NCBI}/1100/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"}
        match = match_gene_annotation_release(candidates, _refseq("GRCh38", REFSEQ_GRCH38_110))
        self.assertEqual([], match.release_ids)

    def test_ambiguous_match_has_no_unique_release(self):
        duplicate_url = f"{NCBI}/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
        candidates = {8: duplicate_url, 9: duplicate_url.replace("https://", "ftp://")}
        match = match_gene_annotation_release(candidates, _refseq("GRCh38", REFSEQ_GRCH38_RS_2025_08))
        self.assertEqual(2, len(match.release_ids))
        self.assertIsNone(match.unique_release_id)
