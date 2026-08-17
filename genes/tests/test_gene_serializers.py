from django.test import TestCase

from annotation.fake_annotation import get_fake_annotation_version
from annotation.tests.test_data_fake_genes import create_fake_transcript_version
from genes.models import TranscriptVersion
from genes.serializers import GeneDetailSerializer, GeneSymbolDetailSerializer
from snpdb.models import GenomeBuild


class TranscriptVersionCoordinatesTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.transcript_version = create_fake_transcript_version(cls.grch37)
        cls.contig = cls.transcript_version.contig

    def test_get_coordinates(self):
        coordinates = self.transcript_version.get_coordinates()
        self.assertIsNotNone(coordinates)
        self.assertEqual(self.contig, coordinates.contig)
        self.assertEqual(self.contig.name, coordinates.chrom)
        self.assertEqual(self.transcript_version.strand, coordinates.strand)
        self.assertIn(coordinates.chrom, coordinates.coordinates)

    def test_get_coordinates_none_for_missing_build(self):
        """ A transcript with no cdot data for a build returns None rather than raising """
        tv = TranscriptVersion(transcript=self.transcript_version.transcript, version=99,
                               gene_version=self.transcript_version.gene_version, genome_build=self.grch37,
                               contig=self.contig, import_source=self.transcript_version.import_source,
                               data={})
        self.assertIsNone(tv.get_coordinates())


class GeneDetailSerializerTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.transcript_version = create_fake_transcript_version(cls.grch37)
        cls.gene = cls.transcript_version.gene_version.gene
        cls.gene_symbol = cls.transcript_version.gene_version.gene_symbol
        cls.contig = cls.transcript_version.contig

    def test_gene_detail_has_versions_and_transcripts(self):
        data = GeneDetailSerializer(self.gene).data
        self.assertEqual(self.gene.identifier, data["identifier"])
        self.assertTrue(data["versions"], "Gene has at least one version")
        version = data["versions"][0]
        self.assertEqual(self.gene_symbol.symbol, version["gene_symbol"])
        contig_ids = {c["id"] for c in version["contigs"]}
        self.assertIn(self.contig.pk, contig_ids)
        transcript_accessions = {tv["accession"] for tv in version["transcript_versions"]}
        self.assertIn(self.transcript_version.accession, transcript_accessions)

    def test_genome_build_context_filters_versions(self):
        grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        data = GeneDetailSerializer(self.gene, context={"genome_build": grch38}).data
        self.assertEqual([], data["versions"], "GRCh38 has no versions for this fake gene")

    def test_transcript_contig_has_id(self):
        """ The contig pk is included so callers can filter variant queries on locus__contig_id """
        data = GeneDetailSerializer(self.gene).data
        transcript = data["versions"][0]["transcript_versions"][0]
        self.assertEqual(self.contig.pk, transcript["contig"]["id"])
        self.assertEqual(self.contig.name, transcript["coordinates"]["chrom"])


class GeneSymbolDetailSerializerTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        get_fake_annotation_version(cls.grch37)
        cls.transcript_version = create_fake_transcript_version(cls.grch37)
        cls.gene_symbol = cls.transcript_version.gene_version.gene_symbol
        cls.contig = cls.transcript_version.contig

    def test_gene_symbol_detail_lists_genes(self):
        data = GeneSymbolDetailSerializer(self.gene_symbol).data
        self.assertEqual(self.gene_symbol.symbol, data["symbol"])
        gene_ids = {g["identifier"] for g in data["genes"]}
        self.assertIn(self.transcript_version.gene_version.gene_id, gene_ids)
