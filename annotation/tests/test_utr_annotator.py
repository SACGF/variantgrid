"""
UTRAnnotator (#2426), added to VariantAnnotationVersions already in use: the VEP command, the migration that
stamps its pin and flags versions for the backfill, and the display parsing.
"""
from importlib import import_module

from django.apps import apps
from django.conf import settings
from django.test import TestCase
from django.test.utils import override_settings

from annotation.enums import VariantAnnotationPipelineType
from annotation.fake_data import (
    get_fake_annotation_settings_dict,
    get_fake_vep_version,
    retire_seeded_annotation_version,
)
from annotation.models import VariantAnnotationVersion
from annotation.models.models_annotation_run import AnnotationRangeLock
from annotation.vep_annotation import get_vep_command
from annotation.vep_plugin_display import utr_annotator_effects
from genes.enums import AnnotationConsortium
from snpdb.models import GenomeBuild
from snpdb.tests.utils.vcf_testing_utils import slowly_create_test_variant

UTR_ANNOTATOR_FILE = "annotation_data/GRCh37/uORF_5UTR_GRCh37_PUBLIC.txt"

_stamp_utr_annotator_and_flag_backfill = import_module(
    "annotation.migrations.0004_utr_annotator")._stamp_utr_annotator_and_flag_backfill


def _settings_with_utr_annotator(vep_version: str = "116") -> dict:
    d = get_fake_annotation_settings_dict(columns_version=5)
    d["ANNOTATION_VEP_VERSION"] = vep_version
    d["ANNOTATION"][settings.BUILD_GRCH37]["vep_config"]["utr_annotator"] = UTR_ANNOTATOR_FILE
    return d


class UTRAnnotatorVepCommandTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")

    def _plugins(self) -> list[str]:
        cmd = get_vep_command("in.vcf", "out.vcf", self.grch37, AnnotationConsortium.ENSEMBL,
                              VariantAnnotationPipelineType.STANDARD)
        return [cmd[i + 1] for i, x in enumerate(cmd) if x == "--plugin"]

    def test_configured_plugin_passes_its_file(self):
        with override_settings(**_settings_with_utr_annotator()):
            plugins = self._plugins()
        self.assertIn(f"UTRAnnotator,file={settings.ANNOTATION_VEP_BASE_DIR}/{UTR_ANNOTATOR_FILE}", plugins)

    def test_unconfigured_plugin_is_skipped(self):
        with override_settings(**get_fake_annotation_settings_dict(columns_version=5)):
            plugins = self._plugins()
        self.assertFalse(any(p.startswith("UTRAnnotator") for p in plugins))


class UTRAnnotatorMigrationTests(TestCase):

    @classmethod
    def setUpTestData(cls):
        cls.grch37 = GenomeBuild.get_name_or_alias("GRCh37")
        retire_seeded_annotation_version(cls.grch37)

    def _vav(self, status, columns_version=5, vep=116, annotated=True) -> VariantAnnotationVersion:
        kwargs = get_fake_vep_version(self.grch37, AnnotationConsortium.ENSEMBL, columns_version)
        kwargs.update(vep=vep, utr_annotator=None, status=status)
        vav = VariantAnnotationVersion.objects.create(**kwargs)
        if annotated:
            variant = slowly_create_test_variant("1", 1000 + vav.pk, "A", "G", self.grch37)
            AnnotationRangeLock.objects.create(version=vav, min_variant=variant, max_variant=variant, count=1)
        return vav

    def _migrate(self, *vavs) -> list[VariantAnnotationVersion]:
        with override_settings(**_settings_with_utr_annotator()):
            _stamp_utr_annotator_and_flag_backfill(apps, None)
        for vav in vavs:
            vav.refresh_from_db()
        return vavs

    def test_active_version_stamped_and_flagged(self):
        vav, = self._migrate(self._vav(VariantAnnotationVersion.Status.ACTIVE))
        self.assertEqual(vav.utr_annotator, "uORF_5UTR_GRCh37_PUBLIC.txt")
        self.assertFalse(vav.backfilled_utr_annotator)

    def test_historical_version_untouched(self):
        vav, = self._migrate(self._vav(VariantAnnotationVersion.Status.HISTORICAL))
        self.assertIsNone(vav.utr_annotator)
        self.assertTrue(vav.backfilled_utr_annotator)

    def test_nothing_to_backfill_stamped_but_not_flagged(self):
        """ A columns_version the plugin's columns don't cover, or a version with no variants yet """
        older_columns = self._vav(VariantAnnotationVersion.Status.ACTIVE, columns_version=4)
        unannotated = self._vav(VariantAnnotationVersion.Status.NEW, annotated=False)
        for vav in self._migrate(older_columns, unannotated):
            self.assertEqual(vav.utr_annotator, "uORF_5UTR_GRCh37_PUBLIC.txt")
            self.assertTrue(vav.backfilled_utr_annotator)


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

