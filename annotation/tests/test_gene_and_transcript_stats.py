from unittest.mock import patch

from django.test import TestCase

from annotation.views import _get_gene_and_transcript_stats
from genes.enums import AnnotationConsortium
from genes.models import CdotDataVersion
from snpdb.models import GenomeBuild

STATS = {"genes": 1, "transcripts": 2, "import_sources": {}}


@patch.object(CdotDataVersion, "count_transcript_stats", return_value=STATS)
class GeneAnnotationStatsTest(TestCase):
    """ The annotation page's gene/transcript counts are stored on CdotDataVersion (#2530) """

    def setUp(self):
        super().setUp()
        self.grch38 = GenomeBuild.grch38()

    def _stats(self):
        return _get_gene_and_transcript_stats(self.grch38, AnnotationConsortium.REFSEQ)

    def test_counted_once_when_missing(self, mock_count):
        cdot_data_version = CdotDataVersion.objects.create(genome_build=self.grch38,
                                                           annotation_consortium=AnnotationConsortium.REFSEQ,
                                                           cdot_version="0.2.34")
        self._stats()
        stats = self._stats()
        self.assertEqual(1, mock_count.call_count)
        self.assertEqual(2, stats["transcripts"])
        self.assertEqual("0.2.34", stats["latest_cdot"])
        cdot_data_version.refresh_from_db()
        self.assertEqual(STATS, cdot_data_version.transcript_stats)
