from django.conf import settings
from django.test import TestCase
from django.test.utils import override_settings

from annotation import vep_columns
from annotation.fake_annotation import get_fake_annotation_settings_dict
from annotation.models import VariantAnnotationVersion
from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.vep_annotation import get_vep_command, get_vep_skipped_variants_filename
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

    def test_skipped_variants_file_named_off_output(self):
        """ #1701: VEP's own list of what it dropped, per-attempt via the output name """
        cmd = self._cmd(AnnotationConsortium.ENSEMBL)
        self.assertIn("--skipped_variants_file", cmd)
        self.assertEqual(cmd[cmd.index("--skipped_variants_file") + 1],
                         get_vep_skipped_variants_filename("out.vcf"))

    @override_settings(ANNOTATION_VEP_MEMORY_LIMIT_GB=4)
    def test_memory_limit_wraps_whole_command(self):
        """ #1710: prlimit has to be outermost so everything VEP starts inherits the heap ceiling """
        cmd = self._cmd(AnnotationConsortium.ENSEMBL)
        self.assertEqual(cmd[:3], ["prlimit", f"--data={4 * 1024 ** 3}", "--"])

    @override_settings(ANNOTATION_VEP_MEMORY_LIMIT_GB=4,
                       ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT="/path/perlbrew_runner.sh")
    def test_memory_limit_wraps_perlbrew_runner(self):
        cmd = self._cmd(AnnotationConsortium.ENSEMBL)
        self.assertEqual(cmd[:4], ["prlimit", f"--data={4 * 1024 ** 3}", "--", "/path/perlbrew_runner.sh"])

    @override_settings(ANNOTATION_VEP_MEMORY_LIMIT_GB=None, ANNOTATION_VEP_PERLBREW_RUNNER_SCRIPT=None)
    def test_no_memory_limit_leaves_command_unwrapped(self):
        cmd = self._cmd(AnnotationConsortium.ENSEMBL)
        self.assertNotIn("prlimit", cmd)
        self.assertTrue(cmd[0].endswith("vep"))

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
    """ columns_version 5 (the helper pins the #1638 plugin data files) at a given VEP version. """
    d = get_fake_annotation_settings_dict(columns_version=5)
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

    def test_open_targets_l2g_scores_grch38_only(self):
        """ #1822 - the raw per-record L2G array is a second column off the same source field """
        for build_name, expected in [("GRCh38", True), ("GRCh37", False)]:
            with self.subTest(build_name):
                columns = vep_columns.visible_columns_for(
                    genome_build_name=build_name, columns_version=5,
                    pipeline_type=VariantAnnotationPipelineType.STANDARD)
                self.assertEqual("open_targets_gwas_l2g_scores" in columns, expected)

    def test_open_targets_l2g_scores_not_in_columns_version_4(self):
        columns = vep_columns.visible_columns_for(
            genome_build_name="GRCh38", columns_version=4,
            pipeline_type=VariantAnnotationPipelineType.STANDARD)
        self.assertNotIn("open_targets_gwas_l2g_scores", columns)

    def test_vep115_omits_eve_and_promoterai(self):
        with override_settings(**_v5_settings(vep_version="115")):
            plugins = self._plugins()
        self.assertTrue(any(p.startswith("ProtVar,db=") for p in plugins))
        self.assertTrue(any(p.startswith("OpenTargets,file=") for p in plugins))
        self.assertFalse(any(p.startswith("EVE,file=") for p in plugins))
        self.assertFalse(any(p.startswith("PromoterAI,file=") for p in plugins))


def _cosmic_settings(cosmic_basename: str) -> dict:
    """ columns_version 5, with the COSMIC custom VCF swapped for a given release. """
    d = get_fake_annotation_settings_dict(columns_version=5)
    d["ANNOTATION"][settings.BUILD_GRCH38]["vep_config"]["cosmic"] = \
        f"annotation_data/GRCh38/{cosmic_basename}"
    return d


class GetVepCommandCosmicTests(TestCase):
    """ #1673 - COSMIC renames the sample count INFO field per release, and VEP silently writes an
        empty CSQ column for a --custom field the file doesn't have, so asking for the wrong name
        leaves cosmic_count null rather than failing. """

    @classmethod
    def setUpTestData(cls):
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")

    def _cosmic_custom_fields(self) -> list[str]:
        cmd = get_vep_command(
            "in.vcf", "out.vcf", self.grch38, AnnotationConsortium.ENSEMBL,
            VariantAnnotationPipelineType.STANDARD,
        )
        customs = [cmd[i + 1] for i, x in enumerate(cmd) if x == "--custom"]
        cosmic = [c for c in customs if "short_name=COSMIC," in c]
        self.assertEqual(len(cosmic), 1)
        params = dict(p.split("=", 1) for p in cosmic[0].split(","))
        return params["fields"].split("%")

    def test_v99_asks_for_sample_count(self):
        with override_settings(**_cosmic_settings("Cosmic_GenomeScreensMutant_v99_GRCh38.vcf.gz")):
            fields = self._cosmic_custom_fields()
        self.assertEqual(sorted(fields), ["LEGACY_ID", "SAMPLE_COUNT"])

    def test_v101_asks_for_genome_screen_sample_count(self):
        with override_settings(**_cosmic_settings("Cosmic_GenomeScreensMutant_Normal_v101_GRCh38.vcf.gz")):
            fields = self._cosmic_custom_fields()
        self.assertEqual(sorted(fields), ["GENOME_SCREEN_SAMPLE_COUNT", "LEGACY_ID"])
