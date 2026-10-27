"""
UTRAnnotator, AVI and SpliceVault (#2426) and the AlphaGenome splicing custom VCF (#2433), added to
VariantAnnotationVersions already in use: the VEP command, enable_vep_plugin_backfill / the migration stamping their
pins and flagging versions for the backfill, the non-commercial data setting, and the display parsing.
"""
import io
import tempfile
from importlib import import_module
from unittest.mock import patch

from django.apps import apps
from django.conf import settings
from django.core.management import CommandError, call_command
from django.db import connection
from django.db.migrations.loader import MigrationLoader
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
from annotation.vep_config import vep_component_version_kwargs
from annotation.vep_plugin_backfill import (
    ENABLE_COMMAND,
    SOURCE_PINS,
    backfill_fork,
    backfill_source_columns,
    backfill_steps,
    backfill_vep_command,
    vep_mismatch_fix,
)
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
ALPHAGENOME_SPLICING_FILE = "annotation_data/GRCh38/alphagenome_splicing_snvs.GRCh38.vcf.gz"

_stamp_pins_and_flag_backfill = import_module(
    "annotation.migrations.0004_utr_annotator_avi_splicevault")._stamp_pins_and_flag_backfill
_backfill_note = import_module(
    "annotation.migrations.0005_utr_annotator_avi_splicevault_manual_backfill")._backfill_note
_flag_versions_annotated_by_old_plugin = import_module(
    "annotation.migrations.0009_utr_effects_backfill")._flag_versions_annotated_by_old_plugin


def _settings_with_plugins(vep_version: str = "116") -> dict:
    d = get_fake_annotation_settings_dict(columns_version=5)
    d["ANNOTATION_VEP_VERSION"] = vep_version
    d["ANNOTATION"][settings.BUILD_GRCH37]["vep_config"]["utr_annotator"] = UTR_ANNOTATOR_FILE
    for build_name, splicevault_file in SPLICEVAULT_FILES.items():
        d["ANNOTATION"][build_name]["vep_config"]["splicevault"] = splicevault_file
        # AVI and AlphaGenome splicing only have GRCh38 data - configured on both to show the build gate drops them
        d["ANNOTATION"][build_name]["vep_config"]["avi"] = AVI_FILE
        d["ANNOTATION"][build_name]["vep_config"]["alphagenome_splicing"] = ALPHAGENOME_SPLICING_FILE
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

    def _customs(self, build_name: str) -> list[str]:
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        cmd = get_vep_command("in.vcf", "out.vcf", genome_build, AnnotationConsortium.REFSEQ,
                              VariantAnnotationPipelineType.STANDARD)
        return [cmd[i + 1] for i, x in enumerate(cmd) if x == "--custom"]

    def test_alphagenome_splicing_custom_matches_alleles_on_grch38(self):
        with override_settings(**_settings_with_plugins()):
            grch37_customs = self._customs("GRCh37")
            grch38_customs = self._customs("GRCh38")
        base_dir = settings.ANNOTATION_VEP_BASE_DIR
        self.assertIn(f"file={base_dir}/{ALPHAGENOME_SPLICING_FILE},type=exact,num_records=all,"
                      f"short_name=AlphaGenome_splicing,fields=score,format=vcf", grch38_customs)
        self.assertFalse(any("AlphaGenome" in c for c in grch37_customs))

    def test_non_commercial_data_left_out_when_not_allowed(self):
        settings_dict = _settings_with_plugins()
        with override_settings(**settings_dict, ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED=False):
            customs = self._customs("GRCh38")
            pins = vep_component_version_kwargs(settings.ANNOTATION["GRCh38"])
        self.assertFalse(any("AlphaGenome" in c for c in customs))
        self.assertIsNone(pins["alphagenome_splicing"])
        self.assertEqual(pins["avi"], "alphagenome_variant_impact_score_snvs.tsv.gz")

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
        kwargs.update(vep=vep, utr_annotator=None, avi=None, splicevault=None, alphagenome_splicing=None,
                      status=status)
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

    def test_earlier_migration_stamps_only_the_pins_its_model_has(self):
        """ 0004's historical VariantAnnotationVersion has no alphagenome_splicing (added by 0007) """
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        state = MigrationLoader(connection).project_state(("annotation", "0004_utr_annotator_avi_splicevault"))
        with override_settings(**_settings_with_plugins()):
            _stamp_pins_and_flag_backfill(state.apps, None)
        vav.refresh_from_db()
        self.assertEqual(vav.avi, "alphagenome_variant_impact_score_snvs.tsv.gz")
        self.assertIsNone(vav.alphagenome_splicing)

    def test_command_refuses_missing_files(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        with tempfile.TemporaryDirectory() as empty_dir, \
                override_settings(**_settings_with_plugins(), ANNOTATION_VEP_BASE_DIR=empty_dir), \
                self.assertRaisesRegex(CommandError, "avi.*MISSING"):
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

    def test_command_stamps_new_and_active_and_flags_only_the_newest(self):
        older = self._vav(self.grch38, VariantAnnotationVersion.Status.NEW)
        historical = self._vav(self.grch38, VariantAnnotationVersion.Status.HISTORICAL)
        active = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        output = self._enable("--genome-build", "GRCh38")
        for vav in (older, historical, active):
            vav.refresh_from_db()
        self.assertEqual(active.splicevault, "SpliceVault_data_GRCh38.refseq_mane.tsv.gz")
        self.assertFalse(active.backfilled_utr_avi_splicevault)
        self.assertIsNone(historical.splicevault)
        self.assertTrue(historical.backfilled_utr_avi_splicevault)
        self.assertEqual(older.avi, "alphagenome_variant_impact_score_snvs.tsv.gz")
        self.assertTrue(older.backfilled_utr_avi_splicevault)
        self.assertIn("Start annotation_workers", output)
        self.assertIn(f"--variant-annotation-version {active.pk} ", output)
        self.assertNotIn(f"--variant-annotation-version {older.pk} ", output)
        self.assertIn("--csq avi_phred=AVI_PHRED", output)
        self.assertEqual(active.alphagenome_splicing, "alphagenome_splicing_snvs.GRCh38.vcf.gz")
        self.assertIn("--csq alphagenome_splicing=AlphaGenome_splicing_score", output)

    def test_unannotated_version_stamped_but_not_flagged(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.NEW, annotated=False)
        self._enable("--genome-build", "GRCh38")
        vav.refresh_from_db()
        self.assertEqual(vav.avi, "alphagenome_variant_impact_score_snvs.tsv.gz")
        self.assertTrue(vav.backfilled_utr_avi_splicevault)

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

    def test_command_does_not_stamp_non_commercial_data_when_not_allowed(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        settings_dict = {**_settings_with_plugins(), "ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED": False}
        output = self._enable("--genome-build", "GRCh38", settings_dict=settings_dict)
        vav.refresh_from_db()
        self.assertIsNone(vav.alphagenome_splicing)
        self.assertIsNotNone(vav.avi)
        self.assertIn("not stamping alphagenome_splicing", output)
        self.assertNotIn("AlphaGenome_splicing_score", output)

    def test_older_columns_version_stamped_but_not_flagged(self):
        vav = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE, columns_version=4)
        self._enable("--genome-build", "GRCh37")
        vav.refresh_from_db()
        self.assertEqual(vav.splicevault, "SpliceVault_data_hg19.refseq_mane.tsv.gz")
        self.assertTrue(vav.backfilled_utr_avi_splicevault)

    def test_mismatch_fix_stamps_pins_only_when_nothing_else_differs(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE, annotated=False)
        with override_settings(**_settings_with_plugins()):
            component_kwargs = vep_component_version_kwargs(settings.ANNOTATION[settings.BUILD_GRCH38])
            pins = {pin_field: component_kwargs[pin_field] for pin_field in SOURCE_PINS}
            self.assertIn(ENABLE_COMMAND, vep_mismatch_fix(vav, pins))
            self.assertIn("create_new_variant_annotation_version", vep_mismatch_fix(vav, {**pins, "vep": 117}))

    def test_backfill_note_imports_only_the_plugins_each_build_covers(self):
        """ AVI is configured on GRCh37 too, but has no GRCh37 columns """
        self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        self._enable("--genome-build", "GRCh37")
        with override_settings(**_settings_with_plugins()):
            note = _backfill_note(apps)
        self.assertIn("--csq splicevault_spliceai_delta=SpliceVault_SpliceAI_delta", note)
        self.assertIn("--csq utr_annotator_consequence=5UTR_consequence", note)
        self.assertNotIn("AVI_PHRED", note)
        self.assertIn(f"--plugin UTRAnnotator,file={settings.ANNOTATION_VEP_BASE_DIR}/{UTR_ANNOTATOR_FILE}", note)
        self.assertIn("enable_vep_plugin_backfill", note)

    def test_backfill_vep_command_keeps_only_the_backfilled_sources(self):
        vav = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        with override_settings(**_settings_with_plugins()):
            self._enable("--genome-build", "GRCh38")
            vav.refresh_from_db()
            columns_by_source = backfill_source_columns(vav)
            cmd = backfill_vep_command(vav, columns_by_source, ["Feature", "PICK", "AVI_PHRED"], "in.vcf.gz",
                                       "out.vcf.gz")
        plugins = [cmd[i + 1].split(",")[0] for i, x in enumerate(cmd) if x == "--plugin"]
        customs = [cmd[i + 1] for i, x in enumerate(cmd) if x == "--custom"]
        self.assertEqual(sorted(plugins), ["AVI", "SpliceVault"])
        self.assertEqual(len(customs), 1)
        self.assertIn("short_name=AlphaGenome_splicing", customs[0])
        for flag in ("--hgvs", "--check_existing", "--sift"):
            self.assertNotIn(flag, cmd)
        self.assertEqual(cmd[-6:], ["--fork", str(backfill_fork()), "--buffer_size", "5000",
                                    "--fields", "Feature,PICK,AVI_PHRED"])
        self.assertEqual(cmd.count("--buffer_size"), 1)

    def _steps(self, vav) -> list[str]:
        with override_settings(**_settings_with_plugins()):
            vav.refresh_from_db()
            return backfill_steps(vav, AnnotationRangeLock)

    def test_utr_effects_migration_flags_versions_that_ran_utr_annotator(self):
        """ #2499: only a NEW/ACTIVE version with the pin and annotated variants ran the old plugin """
        ran = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        unannotated = self._vav(self.grch37, VariantAnnotationVersion.Status.NEW, annotated=False)
        historical = self._vav(self.grch37, VariantAnnotationVersion.Status.HISTORICAL)
        no_utr = self._vav(self.grch38, VariantAnnotationVersion.Status.ACTIVE)
        VariantAnnotationVersion.objects.filter(pk__in=[ran.pk, unannotated.pk, historical.pk]) \
            .update(utr_annotator="uORF_5UTR_GRCh37_PUBLIC.txt")
        _flag_versions_annotated_by_old_plugin(apps, None)
        for vav in (ran, unannotated, historical, no_utr):
            vav.refresh_from_db()
        self.assertFalse(ran.backfilled_utr_effects)
        self.assertTrue(unannotated.backfilled_utr_effects)
        self.assertTrue(historical.backfilled_utr_effects)
        self.assertTrue(no_utr.backfilled_utr_effects)

    def test_utr_effects_pass_alone_reruns_utr_annotator_over_variants_with_a_consequence(self):
        vav = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        self._enable("--genome-build", "GRCh37")
        VariantAnnotationVersion.objects.filter(pk=vav.pk).update(backfilled_utr_avi_splicevault=True,
                                                                    backfilled_utr_effects=False)
        steps = self._steps(vav)
        dump, vep, _dry_run, import_step = steps[1], steps[2], steps[3], steps[4]
        self.assertIn("--columns utr_annotator_consequence,utr_annotator_annotation ", dump)
        self.assertIn("--not-null utr_annotator_consequence", dump)
        self.assertNotIn("--max-variant-id", dump)
        self.assertIn("--plugin UTRAnnotator,", vep)
        self.assertNotIn("SpliceVault", vep)
        self.assertNotIn("splicevault", import_step)
        self.assertTrue(import_step.endswith("--mark-backfilled backfilled_utr_effects"))

    def test_both_pending_is_one_unbounded_pass_marking_both_flags(self):
        vav = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        self._enable("--genome-build", "GRCh37")
        VariantAnnotationVersion.objects.filter(pk=vav.pk).update(backfilled_utr_effects=False)
        steps = self._steps(vav)
        dump, vep, import_step = steps[1], steps[2], steps[4]
        self.assertNotIn("--max-variant-id", dump)
        self.assertNotIn("--not-null", dump)
        self.assertIn("--csq splicevault_spliceai_delta=SpliceVault_SpliceAI_delta", import_step)
        self.assertIn("--plugin UTRAnnotator,", vep)
        self.assertIn("--plugin SpliceVault,", vep)
        self.assertTrue(import_step.endswith(
            "--mark-backfilled backfilled_utr_avi_splicevault,backfilled_utr_effects"))

    def test_plugin_backfill_alone_stays_bounded(self):
        vav = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        self._enable("--genome-build", "GRCh37")
        steps = self._steps(vav)
        self.assertIn("--max-variant-id", steps[1])
        self.assertTrue(steps[4].endswith("--mark-backfilled backfilled_utr_avi_splicevault"))
        VariantAnnotationVersion.objects.filter(pk=vav.pk).update(backfilled_utr_avi_splicevault=True)
        self.assertEqual(self._steps(vav), [])

    def test_mark_backfilled_accepts_several_flags(self):
        vav = self._vav(self.grch37, VariantAnnotationVersion.Status.ACTIVE)
        VariantAnnotationVersion.objects.filter(pk=vav.pk).update(backfilled_utr_avi_splicevault=False,
                                                                    backfilled_utr_effects=False)
        with tempfile.NamedTemporaryFile(suffix=".vcf") as empty_vcf, override_settings(**_settings_with_plugins()):
            empty_vcf.write(b"##fileformat=VCFv4.2\n"
                            b'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Feature|PICK|5UTR_consequence">\n'
                            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            empty_vcf.flush()
            call_command("annotation_backfill_columns", "--import", "--genome-build", "GRCh37",
                         "--variant-annotation-version", str(vav.pk), "--columns", "utr_annotator_consequence",
                         "--csq", "utr_annotator_consequence=5UTR_consequence", "--vcf", empty_vcf.name,
                         "--mark-backfilled", "backfilled_utr_avi_splicevault,backfilled_utr_effects",
                         stdout=io.StringIO())
        vav.refresh_from_db()
        self.assertTrue(vav.backfilled_utr_avi_splicevault)
        self.assertTrue(vav.backfilled_utr_effects)


class UTRAnnotatorEffectsTests(TestCase):

    def test_effects_split_into_records(self):
        """ VEP's VCF output turns the plugin's ',' between records into '&' """
        annotation = ("KozakStrength=Weak:type=uORF:Evidence=NA"
                      "&Evidence=True:Unlisted=1:type=inFrame_oORF:KozakStrength=Strong")
        effects = utr_annotator_effects(annotation)
        self.assertEqual(effects, [
            {"type": "uORF", "KozakStrength": "Weak", "Evidence": "NA"},
            {"type": "inFrame_oORF", "KozakStrength": "Strong", "Evidence": "True", "Unlisted": "1"},
        ])
        # The plugin's key order varies record to record - ours doesn't
        self.assertEqual([list(effect) for effect in effects],
                         [["type", "KozakStrength", "Evidence"], ["type", "KozakStrength", "Evidence", "Unlisted"]])

    def test_no_annotation(self):
        self.assertEqual(utr_annotator_effects(None), [])
