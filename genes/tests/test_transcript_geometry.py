from django.test import TestCase

from annotation.tests.test_data_fake_genes import create_gata2_transcript_version
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import TranscriptGeometry
from snpdb.models.models_genome import GenomeBuild


class TestTranscriptGeometry(TestCase):
    """ GATA2 NM_001145661.2 is minus strand / 7 exons, so the final exon in transcript order
        is the genomically-leftmost one. """

    @classmethod
    def setUpTestData(cls):
        super().setUpTestData()
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.transcript_version = create_gata2_transcript_version(cls.grch37)

    def test_geometry_matches_transcript_version(self):
        geometry = TranscriptGeometry.for_transcript_version(self.transcript_version)
        # start_codon / stop_codon in the cdot blob are cDNA coords, so they pin the derived lengths
        data = self.transcript_version.data
        self.assertEqual(geometry.fivep_utr_length, data["start_codon"])
        self.assertEqual(geometry.coding_length, data["stop_codon"] - data["start_codon"])
        self.assertFalse(geometry.single_exon)

    def test_last_junction_is_cdna_end_of_penultimate_exon(self):
        geometry = TranscriptGeometry.for_transcript_version(self.transcript_version)
        exons = self.transcript_version.genome_build_data["exons"]
        final_exon, penultimate_exon = exons[0], exons[1]  # minus strand: transcript order is reversed
        self.assertEqual(geometry.last_junction_cdna, penultimate_exon[4])

        transcript_length = sum(end - start for start, end, *_ in exons)
        self.assertEqual(geometry.last_junction_cdna,
                         transcript_length - (final_exon[1] - final_exon[0]))

    def test_sparse_cdot_data_returns_none(self):
        self.transcript_version.data["genome_builds"][self.grch37.name].pop("exons")
        self.assertIsNone(TranscriptGeometry.for_transcript_version(self.transcript_version))
