import json
import os
import tempfile

from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from annotation.external_annotation import (
    DUMP_METADATA_SCHEMA_VERSION,
    VariantIdAlignmentError,
    build_dump_metadata,
    dump_external_annotation_runs,
    parse_dump_metadata,
    verify_annotated_vcf_variant_ids,
    write_dump_metadata,
)
from annotation.fake_annotation import (
    create_fake_variants,
    get_fake_annotation_settings_dict,
    get_fake_annotation_version,
)
from annotation.models import AnnotationRangeLock, AnnotationRun
from annotation.models.models import VariantAnnotationVersion
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.tasks.annotate_variants import annotate_variants
from snpdb.models import GenomeBuild, Variant


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class ExternalAnnotationMetadataTests(TestCase):
    """ Step 1 (#1568): self-describing dump filenames + sidecar metadata round-trip. """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.variant_annotation_version = cls.annotation_version.variant_annotation_version
        create_fake_variants(cls.genome_build)

    def _make_run(self):
        variants = list(Variant.objects.order_by("pk"))
        min_variant = variants[0]
        max_variant = variants[-1]
        annotation_range_lock = AnnotationRangeLock.objects.create(
            version=self.variant_annotation_version,
            min_variant=min_variant, max_variant=max_variant, count=len(variants))
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock)
        annotation_run.dump_count = len(variants)
        annotation_run.save()
        return annotation_run, min_variant, max_variant

    def test_dump_filename_is_self_describing(self):
        annotation_run, _, _ = self._make_run()
        basename = os.path.basename(annotation_run.get_dump_filename())
        vav = self.variant_annotation_version
        self.assertIn(self.genome_build.name, basename)
        self.assertIn(f"vep{vav.vep}", basename)
        self.assertIn(f"cv{vav.columns_version}", basename)
        self.assertIn(f"run{annotation_run.pk}", basename)
        self.assertIn("standard", basename)
        self.assertTrue(basename.endswith(".vcf"))

    def test_metadata_filename_is_sidecar(self):
        annotation_run, _, _ = self._make_run()
        dump_filename = annotation_run.get_dump_filename()
        meta_filename = annotation_run.get_dump_metadata_filename()
        self.assertEqual(meta_filename, os.path.splitext(dump_filename)[0] + ".meta.json")

    def test_get_dump_metadata_contents(self):
        annotation_run, min_variant, max_variant = self._make_run()
        meta = build_dump_metadata(annotation_run)
        vav = self.variant_annotation_version

        self.assertEqual(meta["schema"], DUMP_METADATA_SCHEMA_VERSION)
        self.assertEqual(meta["annotation_run_pk"], annotation_run.pk)
        self.assertEqual(meta["pipeline_type"], annotation_run.pipeline_type)
        self.assertEqual(meta["genome_build"], self.genome_build.name)
        self.assertEqual(meta["annotation_consortium"], vav.annotation_consortium)

        vav_meta = meta["variant_annotation_version"]
        self.assertEqual(vav_meta["pk"], vav.pk)
        self.assertEqual(vav_meta["vep"], vav.vep)
        self.assertEqual(vav_meta["columns_version"], vav.columns_version)
        self.assertEqual(vav_meta["gene_annotation_release"]["pk"], vav.gene_annotation_release_id)

        # The range coordinate strings are the Variant string repr (#1568 §4.1 sanity check)
        self.assertEqual(meta["range"]["min_variant"], str(min_variant.coordinate))
        self.assertEqual(meta["range"]["max_variant"], str(max_variant.coordinate))
        self.assertEqual(meta["range"]["count"], annotation_run.annotation_range_lock.count)

        self.assertEqual(meta["batch"]["annotation_vep_batch_min"], 5000)
        self.assertEqual(meta["batch"]["annotation_vep_batch_max"], 25_000)

    def test_metadata_write_parse_round_trip(self):
        annotation_run, _, _ = self._make_run()
        with tempfile.TemporaryDirectory() as tmp_dir:
            with override_settings(ANNOTATION_VCF_DUMP_DIR=tmp_dir):
                meta_filename = write_dump_metadata(annotation_run)
                self.assertTrue(os.path.exists(meta_filename))

                parsed = parse_dump_metadata(meta_filename)
                self.assertEqual(parsed, build_dump_metadata(annotation_run))

    def test_parse_dump_metadata_rejects_bad_schema(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            bad_filename = os.path.join(tmp_dir, "bad.meta.json")
            with open(bad_filename, "w") as f:
                json.dump({"schema": 999}, f)
            with self.assertRaises(ValueError):
                parse_dump_metadata(bad_filename)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class ExternalAnnotationStatusTests(TestCase):
    """ Step 2 (#1568): EXTERNAL_DUMP_COMPLETED state + external runs skipped by annotate_variants. """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.variant_annotation_version = cls.annotation_version.variant_annotation_version
        create_fake_variants(cls.genome_build)

    def _make_external_run(self, dump_count=10):
        variants = list(Variant.objects.order_by("pk"))
        annotation_range_lock = AnnotationRangeLock.objects.create(
            version=self.variant_annotation_version,
            min_variant=variants[0], max_variant=variants[-1], count=len(variants))
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock,
                                                      external=True)
        annotation_run.dump_start = timezone.now()
        annotation_run.dump_end = timezone.now()
        annotation_run.dump_count = dump_count
        annotation_run.save()
        return annotation_run

    def test_external_dumped_run_is_awaiting(self):
        annotation_run = self._make_external_run()
        self.assertEqual(annotation_run.get_status(), AnnotationStatus.EXTERNAL_DUMP_COMPLETED)
        self.assertEqual(annotation_run.status, AnnotationStatus.EXTERNAL_DUMP_COMPLETED)
        # Not a completed state - genuinely waiting on the operator
        self.assertNotIn(annotation_run.status, AnnotationStatus.get_completed_states())

    def test_external_empty_dump_finishes(self):
        # Nothing to annotate - finished even though external
        annotation_run = self._make_external_run(dump_count=0)
        self.assertEqual(annotation_run.status, AnnotationStatus.FINISHED)

    def test_external_run_resumes_once_annotated(self):
        annotation_run = self._make_external_run()
        # --import fills vcf_annotated_filename + the upload path sets upload timestamps
        annotation_run.vcf_annotated_filename = "/tmp/fake.vep_annotated.vcf.gz"
        annotation_run.upload_start = timezone.now()
        annotation_run.upload_end = timezone.now()
        annotation_run.save()
        self.assertEqual(annotation_run.status, AnnotationStatus.FINISHED)

    def test_summary_state(self):
        self.assertEqual(AnnotationStatus.get_summary_state(AnnotationStatus.EXTERNAL_DUMP_COMPLETED),
                         "Awaiting external annotation")

    def test_annotate_variants_skips_external_run(self):
        annotation_run = self._make_external_run()
        # Synchronous celery in tests; must return early without running VEP or erroring
        annotate_variants.apply((annotation_run.pk,)).get()
        annotation_run.refresh_from_db()
        self.assertIsNone(annotation_run.vcf_annotated_filename)
        self.assertIsNone(annotation_run.annotation_start)
        self.assertEqual(annotation_run.status, AnnotationStatus.EXTERNAL_DUMP_COMPLETED)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class ExternalAnnotationDumpTests(TestCase):
    """ Step 3 (#1568): annotation_external --dump creates + dumps + parks runs awaiting annotation. """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.variant_annotation_version = cls.annotation_version.variant_annotation_version
        create_fake_variants(cls.genome_build)
        cls.num_standard_variants = Variant.objects.all().count()

    def test_dump_helper_creates_parks_and_writes_files(self):
        with tempfile.TemporaryDirectory() as output_dir:
            annotation_runs = dump_external_annotation_runs(self.variant_annotation_version, output_dir)
            self.assertTrue(annotation_runs)

            dumped_count = 0
            for annotation_run in annotation_runs:
                annotation_run.refresh_from_db()
                self.assertTrue(annotation_run.external)
                self.assertEqual(annotation_run.status, AnnotationStatus.EXTERNAL_DUMP_COMPLETED)
                # VCF + sidecar metadata written into output_dir
                self.assertTrue(os.path.exists(annotation_run.vcf_dump_filename))
                self.assertEqual(os.path.dirname(annotation_run.vcf_dump_filename), output_dir)
                meta_filename = annotation_run.get_dump_metadata_filename(dump_dir=output_dir)
                self.assertTrue(os.path.exists(meta_filename))
                dumped_count += annotation_run.dump_count

            # Range locks cover every standard variant exactly once
            self.assertEqual(dumped_count, self.num_standard_variants)

    def test_dump_command_requires_new_version(self):
        # fake annotation version is ACTIVE - no NEW version to dump
        with tempfile.TemporaryDirectory() as output_dir:
            with self.assertRaises(CommandError):
                call_command("annotation_external", "--dump", "--genome-build", "GRCh37",
                             "--output-dir", output_dir)

    def test_dump_command_against_new_version(self):
        # Move the fake version to NEW so the command will operate on it
        self.variant_annotation_version.status = VariantAnnotationVersion.Status.NEW
        self.variant_annotation_version.save()

        with tempfile.TemporaryDirectory() as output_dir:
            call_command("annotation_external", "--dump", "--genome-build", "GRCh37",
                         "--output-dir", output_dir)

        annotation_runs = AnnotationRun.objects.filter(
            annotation_range_lock__version=self.variant_annotation_version, external=True)
        self.assertTrue(annotation_runs.exists())
        for annotation_run in annotation_runs:
            self.assertEqual(annotation_run.pipeline_type, VariantAnnotationPipelineType.STANDARD)
            self.assertEqual(annotation_run.status, AnnotationStatus.EXTERNAL_DUMP_COMPLETED)


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class ExternalAnnotationVariantIdCheckTests(TestCase):
    """ Step 4 / §6a (#1568): min/max endpoint variant-id alignment check. """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.variant_annotation_version = cls.annotation_version.variant_annotation_version
        create_fake_variants(cls.genome_build)

    def _make_run(self):
        variants = list(Variant.objects.order_by("pk"))
        annotation_range_lock = AnnotationRangeLock.objects.create(
            version=self.variant_annotation_version,
            min_variant=variants[0], max_variant=variants[-1], count=len(variants))
        annotation_run = AnnotationRun.objects.create(annotation_range_lock=annotation_range_lock)
        annotation_run.dump_count = len(variants)
        annotation_run.save()
        return annotation_run, variants

    def test_passes_on_same_db_dump(self):
        # The endpoints round-trip: meta written from this DB, checked against this DB
        annotation_run, _ = self._make_run()
        meta = build_dump_metadata(annotation_run)
        # Must not raise
        verify_annotated_vcf_variant_ids(annotation_run, meta)

    def test_fails_on_misaligned_max_endpoint(self):
        annotation_run, variants = self._make_run()
        meta = build_dump_metadata(annotation_run)
        # Simulate a diverged/wrong DB: the origin recorded a different coordinate at the max id
        other_coordinate = str(variants[0].coordinate)
        self.assertNotEqual(other_coordinate, meta["range"]["max_variant"])
        meta["range"]["max_variant"] = other_coordinate
        with self.assertRaises(VariantIdAlignmentError):
            verify_annotated_vcf_variant_ids(annotation_run, meta)

    def test_fails_on_misaligned_min_endpoint(self):
        annotation_run, variants = self._make_run()
        meta = build_dump_metadata(annotation_run)
        other_coordinate = str(variants[-1].coordinate)
        self.assertNotEqual(other_coordinate, meta["range"]["min_variant"])
        meta["range"]["min_variant"] = other_coordinate
        with self.assertRaises(VariantIdAlignmentError):
            verify_annotated_vcf_variant_ids(annotation_run, meta)

    def test_fails_on_missing_local_variant(self):
        annotation_run, _ = self._make_run()
        meta = build_dump_metadata(annotation_run)
        # Point the recorded range at a coordinate while the local id maps to nothing: emulate by
        # setting the range lock max id beyond any existing variant via a fresh range lock isn't possible
        # (FK PROTECT), so instead assert the mismatch message references the run + recorded coordinate.
        meta["range"]["max_variant"] = "99:999999 A>T"
        with self.assertRaises(VariantIdAlignmentError) as cm:
            verify_annotated_vcf_variant_ids(annotation_run, meta)
        self.assertIn(str(annotation_run.pk), str(cm.exception))
