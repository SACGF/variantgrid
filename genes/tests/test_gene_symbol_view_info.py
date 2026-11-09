"""The gene symbol page's Summary for a gene with no HGNC record or RefSeq summary: name, location and ideogram from
the gene annotation and its transcripts."""
import copy

from django.contrib.auth.models import User
from django.test import TestCase
from django.urls import reverse

from annotation.fake_data import get_fake_annotation_version
from genes.enums import AnnotationConsortium
from genes.fake_data import (
    create_fake_gene_version,
    create_fake_transcript_version,
    insert_transcript_data,
)
from genes.gene_symbol_view_info import GeneSymbolViewInfo
from snpdb.models.models_genome import GenomeBuild


class GeneSymbolViewInfoTest(TestCase):
    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.user = User.objects.get_or_create(username="gene_symbol_view_info_user")[0]
        cls.grch37 = GenomeBuild.grch37()
        get_fake_annotation_version(cls.grch37)
        cls.transcript_version = create_fake_transcript_version(cls.grch37)  # RUNX1 - no HGNC, no RefSeq summary
        cls.gene_symbol = cls.transcript_version.gene_version.gene_symbol

    def _view_info(self, genome_build=None) -> GeneSymbolViewInfo:
        return GeneSymbolViewInfo(gene_symbol=self.gene_symbol, desired_genome_build=genome_build or self.grch37,
                                  user=self.user)

    def _add_transcript(self, transcript_id: str, contig: str, genome_build=None):
        genome_build = genome_build or self.grch37
        data = copy.deepcopy(self.transcript_version.data)
        data["id"] = transcript_id
        build_data = data["genome_builds"].pop(self.grch37.name)
        build_data["contig"] = contig
        data["genome_builds"] = {genome_build.name: build_data}
        gene_version = create_fake_gene_version(genome_build, "ENSG00000159216", self.gene_symbol.symbol,
                                                AnnotationConsortium.ENSEMBL)
        return insert_transcript_data(genome_build, data, gene_version)

    def test_location_spans_the_transcripts_bands(self):
        """ cdot exons [34787800, 35049344) are 1-based 34787801-35049344, and every band they touch is listed """
        location = self._view_info().gene_location
        self.assertEqual((location.contig.name, location.start, location.end), ("21", 34787801, 35049344))
        first_band, last_band = location.cytobands[0], location.cytobands[-1]
        self.assertLessEqual(first_band.start, location.start)
        self.assertGreaterEqual(last_band.end, location.end)
        self.assertTrue(location.cytoband_label.startswith("21q"))
        self.assertIn(f"band={location.band}", location.contig_url)

    def test_location_on_the_contig_most_transcripts_are_on(self):
        """ A PAR gene's transcripts are on X and Y - the location is the contig with more of them """
        self._add_transcript("ENST90000000001.1", "X")
        self._add_transcript("ENST90000000002.1", "X")
        self.assertEqual(self._view_info().gene_location.contig.name, "X")

    def test_description_prefers_the_page_build(self):
        grch38 = GenomeBuild.grch38()
        gene_versions = [("GRCh37 name", self.transcript_version.gene_version),
                         ("GRCh38 name", self._add_transcript("ENST90000000003.1", "21", grch38).gene_version)]
        for description, gene_version in gene_versions:
            gene_version.description = description
            gene_version.save()
        self.assertEqual(self._view_info(grch38).gene_description, "GRCh38 name")
        self.assertEqual(self._view_info(self.grch37).gene_description, "GRCh37 name")

    def test_summary_card_shown_without_hgnc_or_refseq_summary(self):
        self.client.force_login(self.user)
        url = reverse("view_gene_symbol_genome_build",
                      kwargs={"gene_symbol": self.gene_symbol.symbol, "genome_build_name": self.grch37.name})
        response = self.client.get(url)
        self.assertContains(response, 'id="gene-symbol-ideogram"')
        self.assertContains(response, response.context["gene_location"].cytoband_label)
