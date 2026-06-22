from django.test import TestCase
from django.test.utils import override_settings

from annotation.fake_annotation import get_fake_annotation_settings_dict
from annotation.models import VariantAnnotationVersion
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.vep_annotation import get_vep_command
from genes.models_enums import AnnotationConsortium
from snpdb.models import GenomeBuild


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class GetVepCommandTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

    def _cmd(self, annotation_consortium, vav=None):
        return get_vep_command(
            "in.vcf", "out.vcf", self.grch38, annotation_consortium,
            VariantAnnotationPipelineType.STANDARD,
            variant_annotation_version=vav,
        )

    def test_gencode_primary_ensembl_sets_flag(self):
        vav = VariantAnnotationVersion(
            gencode_subset=VariantAnnotationVersion.GencodeSubset.PRIMARY,
            distance=5000,
        )
        cmd = self._cmd(AnnotationConsortium.ENSEMBL, vav=vav)
        self.assertIn("--gencode_primary", cmd)
        self.assertNotIn("--gencode_basic", cmd)

    def test_gencode_basic_ensembl_sets_flag(self):
        vav = VariantAnnotationVersion(
            gencode_subset=VariantAnnotationVersion.GencodeSubset.BASIC,
            distance=5000,
        )
        cmd = self._cmd(AnnotationConsortium.ENSEMBL, vav=vav)
        self.assertIn("--gencode_basic", cmd)
        self.assertNotIn("--gencode_primary", cmd)

    def test_gencode_primary_refseq_skips_flag(self):
        vav = VariantAnnotationVersion(
            gencode_subset=VariantAnnotationVersion.GencodeSubset.PRIMARY,
            distance=5000,
        )
        cmd = self._cmd(AnnotationConsortium.REFSEQ, vav=vav)
        self.assertNotIn("--gencode_primary", cmd)
        self.assertNotIn("--gencode_basic", cmd)

    def test_gencode_subset_none_skips_flag(self):
        vav = VariantAnnotationVersion(gencode_subset=None, distance=5000)
        cmd = self._cmd(AnnotationConsortium.ENSEMBL, vav=vav)
        self.assertNotIn("--gencode_primary", cmd)
        self.assertNotIn("--gencode_basic", cmd)

    def test_distance_read_from_vav_ignores_setting(self):
        vav = VariantAnnotationVersion(distance=10000, gencode_subset=None)
        with override_settings(ANNOTATION_VEP_DISTANCE=5000):
            cmd = self._cmd(AnnotationConsortium.ENSEMBL, vav=vav)
        self.assertIn("--distance", cmd)
        self.assertEqual(cmd[cmd.index("--distance") + 1], "10000")

    def test_no_vav_omits_distance_flag(self):
        with override_settings(ANNOTATION_VEP_DISTANCE=5000):
            cmd = self._cmd(AnnotationConsortium.ENSEMBL, vav=None)
        self.assertNotIn("--distance", cmd)
        self.assertNotIn("--gencode_primary", cmd)
        self.assertNotIn("--gencode_basic", cmd)
