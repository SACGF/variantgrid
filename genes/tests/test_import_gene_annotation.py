import gzip
import io
import json

from django.test import TestCase

from genes.enums import AnnotationConsortium, GeneSymbolAliasSource
from genes.management.commands.import_gene_annotation import Command as GeneAnnotationCommand
from genes.models import (
    CdotDataVersion,
    Gene,
    GeneSymbol,
    GeneSymbolAlias,
    TranscriptVersion,
)
from snpdb.models.models_genome import GenomeBuild


class ImportCdotDataTest(TestCase):
    URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.25_GRCh37.p13_genomic.gff.gz"

    def setUp(self):
        self.genome_build = GenomeBuild.grch37()

    def _cdot_file(self, cdot_version: str, gene_symbol="RUNX1", exons=None) -> io.BytesIO:
        exons = exons or [[36164432, 36164908, 0, 1, 477, None]]

        def _transcript(gene_version):
            return {
                "gene_name": gene_symbol,
                "gene_version": gene_version,
                "biotype": ["mRNA"],
                "genome_builds": {
                    self.genome_build.name: {"url": self.URL, "contig": "NC_000021.8", "strand": "-",
                                             "exons": exons}
                },
            }

        cdot_data = {
            "cdot_version": cdot_version,
            "genes": {
                "861": {"gene_symbol": "RUNX1", "url": self.URL, "biotype": ["protein_coding"]},
                "2623": {"gene_symbol": "GATA2", "url": self.URL, "biotype": ["protein_coding"]},
            },
            "transcripts": {
                "NM_001754.4": _transcript("861"),
                "NM_032638.4": _transcript("2623") | {"gene_name": "GATA2"},
            },
        }
        return io.BytesIO(gzip.compress(json.dumps(cdot_data).encode()))

    def _import(self, cdot_file):
        with gzip.open(cdot_file, "rb") as f:
            GeneAnnotationCommand.import_cdot_data_file(self.genome_build, AnnotationConsortium.REFSEQ, f,
                                                        GeneAnnotationCommand.read_cdot_version(f))

    def _cdot_by_accession(self) -> dict[str, str]:
        return {tv.accession: tv.modified_cdot_version for tv in TranscriptVersion.objects.filter(genome_build=self.genome_build)}

    def test_only_changed_transcripts_restamped(self):
        self._import(self._cdot_file("0.2.34"))
        self._import(self._cdot_file("0.2.35"))
        self.assertEqual({"NM_001754.4": "0.2.34", "NM_032638.4": "0.2.34"}, self._cdot_by_accession())
        self.assertEqual("0.2.35", CdotDataVersion.get_version(self.genome_build, AnnotationConsortium.REFSEQ))

        # Only GATA2 transcript changes
        cdot_data = json.loads(gzip.decompress(self._cdot_file("0.2.36").getvalue()))
        gata2_build_data = cdot_data["transcripts"]["NM_032638.4"]["genome_builds"][self.genome_build.name]
        gata2_build_data["exons"] = [[36164432, 36164999, 0, 1, 568, None]]
        self._import(io.BytesIO(gzip.compress(json.dumps(cdot_data).encode())))

        self.assertEqual({"NM_001754.4": "0.2.34", "NM_032638.4": "0.2.36"}, self._cdot_by_accession())
        gata2_tv = TranscriptVersion.objects.get(genome_build=self.genome_build, transcript_id="NM_032638", version=4)
        self.assertEqual(gata2_build_data["exons"], gata2_tv.genome_build_data["exons"])
        self.assertEqual("0.2.36", CdotDataVersion.get_version(self.genome_build, AnnotationConsortium.REFSEQ))

    def test_transcript_stats_counted(self):
        self._import(self._cdot_file("0.2.34"))
        cdot_data_version = CdotDataVersion.objects.get(genome_build=self.genome_build,
                                                        annotation_consortium=AnnotationConsortium.REFSEQ)
        self.assertEqual({"genes": 2, "transcripts": 2, "import_sources": {self.URL: 2}},
                         cdot_data_version.transcript_stats)

    def test_mitochondrial_fake_transcript_imported(self):
        """ cdot's versionless 'fake-rna-ND4' is stored as version 1 so MT annotation can link to it, and
            matched on re-import (#2139) """
        cdot_data = json.loads(gzip.decompress(self._cdot_file("0.2.34").getvalue()))
        cdot_data["genes"]["4538"] = {"gene_symbol": "ND4", "url": self.URL, "biotype": ["protein_coding"]}
        cdot_data["transcripts"]["fake-rna-ND4"] = {
            "gene_name": "ND4",
            "gene_version": "4538",
            "biotype": ["mRNA"],
            "genome_builds": {
                self.genome_build.name: {"url": self.URL, "contig": "NC_012920.1", "strand": "+",
                                         "exons": [[10759, 12137, 0, 1, 1378, None]]}
            },
        }
        cdot_file = gzip.compress(json.dumps(cdot_data).encode())
        self._import(io.BytesIO(cdot_file))
        self._import(io.BytesIO(cdot_file))  # Re-import matches the existing fake transcript version
        tv = TranscriptVersion.objects.get(genome_build=self.genome_build, transcript_id="fake-rna-ND4")
        self.assertEqual(1, tv.version)
        self.assertTrue(tv.is_cdot_fake)

    def _ncbi_aliases(self) -> set[tuple[str, str]]:
        return set(GeneSymbolAlias.objects.filter(source=GeneSymbolAliasSource.NCBI).values_list("alias", "gene_symbol_id"))

    def test_gene_summaries_and_ncbi_aliases(self):
        """ Summaries update on existing genes; NCBI aliases skip gene symbols and HGNC aliases, and one cdot drops
            is kept """
        GeneSymbol.objects.create(symbol="GATA2")
        GeneSymbolAlias.objects.create(alias="CBFA2", gene_symbol_id="GATA2", source=GeneSymbolAliasSource.HGNC)
        cdot_data = json.loads(gzip.decompress(self._cdot_file("0.2.34").getvalue()))
        cdot_data["genes"]["861"] |= {"summary": "old", "aliases": "AML1, CBFA2, GATA2, PEBP2aB"}
        self._import(io.BytesIO(gzip.compress(json.dumps(cdot_data).encode())))
        self.assertEqual({("AML1", "RUNX1"), ("PEBP2aB", "RUNX1")}, self._ncbi_aliases())

        cdot_data["genes"]["861"] |= {"summary": "new", "aliases": "PEBP2AB, CBF2alpha"}
        self._import(io.BytesIO(gzip.compress(json.dumps(cdot_data).encode())))
        expected = {("AML1", "RUNX1"), ("PEBP2aB", "RUNX1"), ("CBF2alpha", "RUNX1")}
        self.assertEqual(expected, self._ncbi_aliases())
        self.assertEqual("new", Gene.objects.get(pk="861").summary)
