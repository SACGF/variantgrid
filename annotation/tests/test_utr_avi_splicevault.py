"""
UTRAnnotator, AVI and SpliceVault (#2426), added to VariantAnnotationVersions already in use: the VEP command,
enable_vep_plugin_backfill / the migration stamping their pins and flagging versions for the backfill, and the
display parsing.
"""
import io
from importlib import import_module
from unittest.mock import patch

from django.apps import apps
from django.conf import settings
from django.core.management import CommandError, call_command
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.enums import VariantAnnotationPipelineType
from annotation.fake_data import (
    get_fake_annotation_settings_dict,
    get_fake_vep_version,
    retire_seeded_annotation_version,
)
from annotation.models import VariantAnnotationVersion
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.vep_annotation import get_vep_command
from annotation.vep_plugin_display import utr_annotator_effects
from genes.enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

UTR_ANNOTATOR_FILE = "annotation_data/GRCh37/uORF_5UTR_GRCh37_PUBLIC.txt"
SPLICEVAULT_FILES = {
    settings.BUILD_GRCH37: "annotation_data/GRCh37/SpliceVault_data_hg19.refseq_mane.tsv.gz",
    settings.BUILD_GRCH38: "annotation_data/GRCh38/SpliceVault_data_GRCh38.refseq_mane.tsv.gz",
}
AVI_FILE = "annotation_data/GRCh38/alphagenome_variant_impact_score_snvs.tsv.gz"

_stamp_pins_and_flag_backfill = import_module(
    "annotation.migrations.0004_utr_annotator_avi_splicevault")._stamp_pins_and_flag_backfill
_backfill_note = import_module(
    "annotation.migrations.0005_utr_annotator_avi_splicevault_manual_backfill")._backfill_note


def _settings_with_plugins(vep_version: str = "116") -> dict:
    d = get_fake_annotation_settings_dict(columns_version=5)
    d["ANNOTATION_VEP_VERSION"] = vep_version
    d["ANNOTATION"][settings.BUILD_GRCH37]["vep_config"]["utr_annotator"] = UTR_ANNOTATOR_FILE
    for build_name, splicevault_file in SPLICEVAULT_FILES.items():
        d["ANNOTATION"][build_name]["vep_config"]["splicevault"] = splicevault_file
        # AVI only has GRCh38 data - configured on both to show the build gate drops it
        d["ANNOTATION"][build_name]["vep_config"]["avi"] = AVI_FILE
    return d


class PluginVepCommandTests(TestCase):

    def _plugins(self, build_name: str) -> list[str]:
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        cmd = get_vep_command("in.vcf", "out.vcf", genome_build, AnnotationConsortium.REFSEQ,
                              VariantAnnotationPipelineType.STANDARD)
        return [cmd[i + 1] for i, x in enumerate(cmd) if x == "--plugin"]

    def test_configured_plugins_pass_their_files(self):
        with override_settings(**_settings_with_plugins()):
            grch37_plugins = self._plugins("GRCh37")
            grch38_plugins = self._plugins("GRCh38")
        base_dir = settings.ANNOTATION_VEP_BASE_DIR
        self.assertIn(f"UTRAnnotator,file={base_dir}/{UTR_ANNOTATOR_FILE}", grch37_plugins)
        self.assertIn(f"SpliceVault,file={base_dir}/{SPLICEVAULT_FILES['GRCh37']}", grch37_plugins)
        self.assertFalse(any(p.startswith("AVI") for p in grch37_plugins))
        self.assertIn(f"AVI,file={base_dir}/{AVI_FILE},raw=1", grch38_plugins)

    def test_avi_needs_vep_116(self):
        with override_settings(**_settings_with_plugins(vep_version="115")):
            plugins = self._plugins("GRCh38")
        self.assertFalse(any(p.startswith("AVI") for p in plugins))
        self.assertTrue(any(p.startswith("SpliceVault") for p in plugins))

    def test_unconfigured_plugins_are_skipped(self):
        with override_settings(**get_fake_annotation_settings_dict(columns_version=5)):
            plugins = self._plugins("GRCh38")
        self.assertFalse(any(p.startswith(("UTRAnnotator", "AVI", "SpliceVault")) for p in plugins))


class PluginPinStampTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        cls.grch38 = GenomeBuild.get_name_or_alias("GRCh38")
        retire_seeded_annotation_version(cls.grch37)
        retire_seeded_annotation_version(cls.grch38)

    def _vav(self, genome_build, status, columns_version=5, vep=116, annotated=True) -> VariantAnnotationVersion:
        kwargs = get_fake_vep_version(genome_build, AnnotationConsortium.REFSEQ, columns_version)
        kwargs.update(vep=vep, utr_annotator=None, avi=None, splicevault=None, status=status)
        vav = VariantAnnotationVersion.objects.create(**kwargs)
        if annotated:
            variant = slowly_create_test_variant("1", 1000 + vav.pk, "A", "G", genome_build)
            AnnotationRangeLock.objects.create(version=vav, min_variant=variant, max_variant=variant, count=1)
        return vav

    def _enable(self, *args, settings_dict=None) -> str:
        """ The command, with the data files taken as installed """
        stdout = io.StringIO()
        with override_settings(**(settings_dict or _settings_with_plugins())), \
                patch("annotation.management.commands.enable_vep_plugin_backfill.get_missing_vep_files",
                      return_value={}):
            call_command("enable_vep_plugin_backfill", *args, stdout=stdout)
        return stdout.getvalue()

    def test_migration_stamps_from_settings_set_before_migrating(self):
        grch37 = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        grch38 = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        with override_settings(**_settings_with_plugins()):
            _stamp_pins_and_flag_backfill(apps, None)
        grch37.refresh_from_db()
        grch38.refresh_from_db()
        self.assertEqual(grch37.utr_annotator, "uORF_5UTR_GRCh37_PUBLIC.txt")
        self.assertEqual(grch37.splicevault, "SpliceVault_data_hg19.refseq_mane.tsv.gz")
        self.assertEqual(grch38.avi, "alphagenome_variant_impact_score_snvs.tsv.gz")
        self.assertFalse(grch37.backfilled_utr_avi_splicevault)
        self.assertFalse(grch38.backfilled_utr_avi_splicevault)

    def test_command_refuses_missing_files(self):
        """ The test settings' data files aren't installed """
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        with override_settings(**_settings_with_plugins()), self.assertRaisesRegex(CommandError, "avi.*MISSING"):
            call_command("enable_vep_plugin_backfill", "--genome-build", "GRCh38", stdout=io.StringIO())
        vav.refresh_from_db()
        self.assertIsNone(vav.avi)

    def test_command_refuses_runs_vep_annotated_but_not_imported(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        AnnotationRun.objects.create(annotation_range_lock=vav.annotationrangelock_set.get(),
                                     annotation_start=timezone.now(), annotation_end=timezone.now())
        with self.assertRaisesRegex(CommandError, "1 STANDARD annotation runs"):
            self._enable("--genome-build", "GRCh38")
        self._enable("--genome-build", "GRCh38", "--force")
        vav.refresh_from_db()
        self.assertIsNotNone(vav.avi)

    def test_command_stamps_and_flags_new_and_active_only(self):
        active = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        historical = self._vav(self.grch38, VariantAnnotationVersion.Status.HISTORICAL)
        unannotated = self._vav(self.grch38, VariantAnnotationVersion.Status.NEW, annotated=False)
        output = self._enable("--genome-build", "GRCh38")
        for vav in (active, historical, unannotated):
            vav.refresh_from_db()
        self.assertEqual(active.splicevault, "SpliceVault_data_GRCh38.refseq_mane.tsv.gz")
        self.assertFalse(active.backfilled_utr_avi_splicevault)
        self.assertIsNone(historical.splicevault)
        self.assertTrue(historical.backfilled_utr_avi_splicevault)
        self.assertEqual(unannotated.avi, "alphagenome_variant_impact_score_snvs.tsv.gz")
        self.assertTrue(unannotated.backfilled_utr_avi_splicevault)  # No variants to backfill
        self.assertIn("Restart annotation_workers", output)
        self.assertIn(f"--variant-annotation-version {active.pk} ", output)
        self.assertIn("--csq avi_phred=AVI_PHRED", output)

    def test_command_rerun_after_backfill_changes_nothing(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        self._enable("--genome-build", "GRCh38")
        VariantAnnotationVersion.objects.filter(pk=vav.pk).update(backfilled_utr_avi_splicevault=True)
        output = self._enable("--genome-build", "GRCh38")
        vav.refresh_from_db()
        self.assertTrue(vav.backfilled_utr_avi_splicevault)
        self.assertIn("already matches settings", output)
        self.assertNotIn("Restart", output)

    def test_command_dry_run_changes_nothing(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        output = self._enable("--genome-build", "GRCh38", "--dry-run")
        vav.refresh_from_db()
        self.assertIsNone(vav.avi)
        self.assertTrue(vav.backfilled_utr_avi_splicevault)
        self.assertIn("Would stamp", output)
        self.assertIn("avi: None -> alphagenome_variant_impact_score_snvs.tsv.gz", output)

    def test_one_plugin_is_enough_to_flag(self):
        """ GRCh38 with only AVI configured """
        settings_dict = get_fake_annotation_settings_dict(columns_version=5)
        settings_dict["ANNOTATION"][settings.BUILD_GRCH38]["vep_config"]["avi"] = AVI_FILE
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        self._enable("--genome-build", "GRCh38", settings_dict=settings_dict)
        vav.refresh_from_db()
        self.assertIsNone(vav.splicevault)
        self.assertFalse(vav.backfilled_utr_avi_splicevault)

    def test_older_columns_version_stamped_but_not_flagged(self):
        vav = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE, columns_version=4)
        self._enable("--genome-build", "GRCh37")
        vav.refresh_from_db()
        self.assertEqual(vav.splicevault, "SpliceVault_data_hg19.refseq_mane.tsv.gz")
        self.assertTrue(vav.backfilled_utr_avi_splicevault)

    def test_backfill_note_imports_only_the_plugins_each_build_covers(self):
        """ AVI is configured on GRCh37 too, but has no GRCh37 columns """
        self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        self._enable("--genome-build", "GRCh37")
        note = _backfill_note(apps)
        self.assertIn("--csq splicevault_spliceai_delta=SpliceVault_SpliceAI_delta", note)
        self.assertIn("--csq utr_annotator_consequence=5UTR_consequence", note)
        self.assertNotIn("AVI_PHRED", note)
        self.assertIn("except UTRAnnotator, SpliceVault removed", note)
        self.assertIn("enable_vep_plugin_backfill", note)


class UTRAnnotatorEffectsTests(TestCase):

    def test_effects_split_into_records(self):
        """ VEP's VCF output turns the plugin's ',' between records into '&' """
        annotation = "KozakStrength=Weak:type=uORF:Evidence=NA&KozakStrength=Strong:type=inFrame_oORF:Evidence=True"
        self.assertEqual(utr_annotator_effects(annotation), [
            {"KozakStrength": "Weak", "type": "uORF", "Evidence": "NA"},
            {"KozakStrength": "Strong", "type": "inFrame_oORF", "Evidence": "True"},
        ])

    def test_no_annotation(self):
        self.assertEqual(utr_annotator_effects(None), [])
