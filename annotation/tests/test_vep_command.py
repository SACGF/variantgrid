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


def _v5_settings(vep_version: str) -> dict:
    """ columns_version 5 with the #1638 plugin data files configured on GRCh38. """
    d = get_fake_annotation_settings_dict(columns_version=5)
    grch38_cfg = d["ANNOTATION"]["GRCh38"]["vep_config"]
    grch38_cfg.update({
        "protvar": "annotation_data/all_builds/ProtVar_data.db",
        "open_targets": "annotation_data/GRCh38/open_targets_26.03_vep.tsv.bgz",
        "eve": "annotation_data/GRCh38/eve_merged.vcf.gz",
        "popeve": "annotation_data/GRCh38/grch38_popEVE_ukbb_20250715.vcf.gz",
        "promoter_ai": "annotation_data/GRCh38/promoterAI_tss500.tsv.bgz",
    })
    d["ANNOTATION_VEP_VERSION"] = vep_version
    return d


class GetVepCommandColumnsVersion5Tests(TestCase):
    """ #1638 - ProtVar / OpenTargets / EVE / PromoterAI plugin wiring at columns_version 5. """

    @classmethod
    def setUpTestData(cls):
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

    def _plugins(self):
        cmd = get_vep_command(
            "in.vcf", "out.vcf", self.grch38, AnnotationConsortium.ENSEMBL,
            VariantAnnotationPipelineType.STANDARD,
        )
        return [cmd[i + 1] for i, x in enumerate(cmd) if x == "--plugin"]

    def test_vep116_includes_all_four_plugins(self):
        with override_settings(**_v5_settings(vep_version="116")):
            plugins = self._plugins()
        self.assertTrue(any(p.startswith("ProtVar,db=") for p in plugins))
        self.assertTrue(any(p.startswith("OpenTargets,file=") and "cols=all" in p for p in plugins))
        self.assertTrue(any(p.startswith("EVE,file=") and "popeve_file=" in p for p in plugins))
        self.assertTrue(any(p.startswith("PromoterAI,file=") for p in plugins))

    def test_vep115_omits_eve_and_promoterai(self):
        with override_settings(**_v5_settings(vep_version="115")):
            plugins = self._plugins()
        self.assertTrue(any(p.startswith("ProtVar,db=") for p in plugins))
        self.assertTrue(any(p.startswith("OpenTargets,file=") for p in plugins))
        self.assertFalse(any(p.startswith("EVE,file=") for p in plugins))
        self.assertFalse(any(p.startswith("PromoterAI,file=") for p in plugins))
