"""
Classification.get_resolved_hgvs for a build without a stored value (T2T), generated on the fly
"""
from unittest.mock import patch

from django.test import TestCase

from classification.models.classification import Classification
from genes.fake_data import create_gene_level_variant
from snpdb.gene_level_variants import GENE_LEVEL_CONTIG_NAME, GENE_LEVEL_REF, GENE_LEVEL_SVLEN
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_variant import VariantCoordinate


class ClassificationResolvedHGVSTest(TestCase):

    def test_gene_level_variant_has_no_resolved_hgvs_and_reports_nothing(self):
        variant = create_gene_level_variant(VariantCoordinate(chrom=GENE_LEVEL_CONTIG_NAME, position=3236,
                                                              ref=GENE_LEVEL_REF, alt="<GAIN:HGNC:3236>",
                                                              svlen=GENE_LEVEL_SVLEN))
        classification = Classification(evidence={})
        with patch.object(Classification, "get_variant_for_build", return_value=variant), \
                patch("classification.models.classification.report_exc_info") as mock_report_exc_info:
            self.assertIsNone(classification.get_resolved_hgvs(GenomeBuild.t2tv2()))
        mock_report_exc_info.assert_not_called()
