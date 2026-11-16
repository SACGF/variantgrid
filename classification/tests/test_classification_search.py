from unittest.mock import patch

from django.test import SimpleTestCase, TestCase, override_settings

from annotation.models import VariantAnnotation
from classification.enums import ShareLevel, SpecialEKeys, SubmissionSource
from classification.models.classification import Classification
from classification.models.classification_variant_info_models import ImportedAlleleInfo
from classification.signals.classification_search import (
    _variant_location_summary,
    variant_preview_classifications_extra,
)
from classification.tests.models.test_utils import ClassificationTestUtils
from library.preview_request import PreviewKeyValue
from snpdb.enums import AlleleConversionTool, AlleleOrigin
from snpdb.models import Allele, GenomeBuild, GenomeBuildPatchVersion, VariantAllele
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant


class VariantLocationSummaryTestCase(SimpleTestCase):
    """ Where an unclassified variant with no canonical c.HGVS lies, shown in its search preview """

    def test_in_gene_without_transcript(self):
        """ RefSeq annotation of MT genes has a symbol but no transcript """
        va = VariantAnnotation(consequence="frameshift_variant", symbol="ND4")
        self.assertEqual(_variant_location_summary(va), "in ND4")

    def test_upstream(self):
        va = VariantAnnotation(consequence="upstream_gene_variant", distance=3709, symbol="CYTB")
        self.assertEqual(_variant_location_summary(va), "3709bp upstream of CYTB")

    def test_intergenic(self):
        va = VariantAnnotation(consequence="intergenic_variant")
        self.assertEqual(_variant_location_summary(va), "intergenic")


@override_settings(CLASSIFICATION_ALLOW_UNKNOWN_KEYS=True)
class VariantPreviewClassificationsExtraTest(TestCase):
    """ The search preview of a variant: classification counts by allele, c.HGVS for the variant's build """

    def setUp(self):
        ClassificationTestUtils.setUp()
        self.lab, self.user = ClassificationTestUtils.lab_and_user()

    def test_build_without_cached_resolved_hgvs(self):
        """ T2T has no cached resolved c.HGVS column (#2821), the allele's classifications are still counted """
        t2t = GenomeBuild.t2tv2()
        variant = slowly_create_test_variant("1", 1000, "A", "T", t2t)
        allele = Allele.objects.create()
        VariantAllele.objects.create(variant=variant, genome_build=t2t, allele=allele,
                                     origin=AlleleOrigin.IMPORTED_TO_DATABASE,
                                     allele_linking_tool=AlleleConversionTool.SAME_CONTIG)
        data = {SpecialEKeys.GENE_SYMBOL: {"value": "RUNX1"}, SpecialEKeys.ALLELE_ORIGIN: {"value": "germline"}}
        classification = Classification.create(user=self.user, lab=self.lab, data=data, save=True,
                                               source=SubmissionSource.API)
        classification.publish_latest(user=self.user, share_level=ShareLevel.ALL_USERS)
        allele_info = ImportedAlleleInfo.objects.create(
            imported_genome_build_patch_version=GenomeBuildPatchVersion.get_unspecified_patch_version_for(t2t),
            imported_c_hgvs="NM_001754.5:c.1A>T", allele=allele)
        Classification.objects.filter(pk=classification.pk).update(allele_info=allele_info)

        with patch("classification.signals.classification_search._variant_hgvs_extra",
                   return_value=[PreviewKeyValue(None, "variant HGVS", dedicated_row=True)]) as mock_variant_hgvs_extra:
            extras = variant_preview_classifications_extra(None, user=self.user, obj=variant)
        mock_variant_hgvs_extra.assert_called_once_with(variant, t2t)
        self.assertEqual([(extra.key, extra.value) for extra in extras],
                         [("Germline Classification Count", 1), (None, "variant HGVS")])
