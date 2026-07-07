import gzip
import json
import os
import tempfile

from django.conf import settings
from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase
from django.test.utils import override_settings
from django.utils import timezone

from unittest.mock import patch

from annotation.external_annotation import (
    ANNOTATED_VCF_SUFFIX,
    DUMP_METADATA_SCHEMA_VERSION,
    VariantIdAlignmentError,
    build_dump_metadata,
    build_snakemake_config,
    build_vep_command_template,
    dump_external_annotation_runs,
    find_annotated_vcf,
    find_matching_annotation_run,
    find_matching_variant_annotation_version,
    import_external_annotation_runs,
    parse_dump_metadata,
    verify_annotated_vcf_variant_ids,
    write_dump_metadata,
    write_snakemake_bundle,
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

    def test_external_dump_filename_is_gzipped(self):
        annotation_run, _, _ = self._make_run()
        annotation_run.external = True
        annotation_run.save()
        dump_filename = annotation_run.get_dump_filename()
        self.assertTrue(dump_filename.endswith(".vcf.gz"))
        # Sidecar shares the stem (no .vcf/.gz) so the Snakefile derives sample -> {stem}.vcf.gz
        meta_filename = annotation_run.get_dump_metadata_filename()
        self.assertEqual(dump_filename, meta_filename[:-len(".meta.json")] + ".vcf.gz")

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
        # Short label for the status-graph legend; the full status label stays "Awaiting external annotation"
        self.assertEqual(AnnotationStatus.get_summary_state(AnnotationStatus.EXTERNAL_DUMP_COMPLETED),
                         "External")

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
                # VCF + sidecar metadata written into output_dir - external dumps are gzipped
                self.assertTrue(os.path.exists(annotation_run.vcf_dump_filename))
                self.assertTrue(annotation_run.vcf_dump_filename.endswith(".vcf.gz"))
                with gzip.open(annotation_run.vcf_dump_filename, "rt") as f:
                    self.assertTrue(f.readline().startswith("#"))  # valid gzip + VCF header
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


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class ExternalAnnotationSnakemakeTests(TestCase):
    """ Step 6 (#1568 §4.3): the --dump command emits a self-contained Snakemake bundle. """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.variant_annotation_version = cls.annotation_version.variant_annotation_version
        create_fake_variants(cls.genome_build)

    def test_command_template_uses_config_placeholders(self):
        template = build_vep_command_template(self.variant_annotation_version)
        joined = " ".join(template)
        # input/output are wildcards; server paths are rewritten to config placeholders
        self.assertIn("{input}", template)
        self.assertIn("{output}", template)
        self.assertIn("{annotation_base_dir}", joined)
        self.assertIn("{vep_code_dir}", joined)
        # fork + buffer_size are overridable from config even if the server ran fork=1
        self.assertIn("{vep_fork}", joined)
        self.assertIn("{vep_buffer_size}", joined)
        # No raw server path leaks through (everything under the base dir became a placeholder)
        self.assertNotIn(settings.ANNOTATION_BASE_DIR, joined)

    def test_write_snakemake_bundle(self):
        with tempfile.TemporaryDirectory() as output_dir:
            snakefile_path, config_path = write_snakemake_bundle(output_dir, self.variant_annotation_version)
            self.assertTrue(os.path.exists(snakefile_path))
            self.assertTrue(os.path.exists(config_path))

            with open(snakefile_path) as f:
                snakefile = f.read()
            self.assertIn('configfile: "config.yaml"', snakefile)
            self.assertIn("rule vep:", snakefile)
            self.assertIn("rule all:", snakefile)
            self.assertIn(ANNOTATED_VCF_SUFFIX, snakefile)
            # The embedded template must be valid JSON referencing input/output
            self.assertIn("{input}", snakefile)
            self.assertIn("{output}", snakefile)

            with open(config_path) as f:
                config_yaml = f.read()
            for key in ("dump_dir", "output_dir", "annotation_base_dir", "vep_code_dir",
                        "vep_fork", "vep_buffer_size"):
                self.assertIn(key, config_yaml)

    def test_dump_helper_emits_bundle(self):
        with tempfile.TemporaryDirectory() as output_dir:
            dump_external_annotation_runs(self.variant_annotation_version, output_dir)
            self.assertTrue(os.path.exists(os.path.join(output_dir, "Snakefile")))
            self.assertTrue(os.path.exists(os.path.join(output_dir, "config.yaml")))

    def test_template_renders_against_config(self):
        # Simulate what the Snakefile does at runtime: token.format(**config, input=, output=). Every
        # placeholder must resolve, leaving no stray braces (which would crash str.format on the compute box).
        template = build_vep_command_template(self.variant_annotation_version)
        config = build_snakemake_config()
        fmt = dict(config)
        fmt["input"] = "in.vcf"
        fmt["output"] = "out.vcf.gz"
        rendered = [token.format(**fmt) for token in template]
        self.assertIn("in.vcf", rendered)
        self.assertIn("out.vcf.gz", rendered)
        leftover = [token for token in rendered if "{" in token or "}" in token]
        self.assertEqual(leftover, [])


@override_settings(**get_fake_annotation_settings_dict(columns_version=2))
class ExternalAnnotationImportTests(TestCase):
    """ Step 5 (#1568 §4.4): match annotated VCFs back to local runs, §6a pre-flight, upload-only import. """

    @classmethod
    def setUpTestData(cls):
        cls.genome_build = GenomeBuild.get_name_or_alias("GRCh37")
        cls.annotation_version = get_fake_annotation_version(cls.genome_build)
        cls.variant_annotation_version = cls.annotation_version.variant_annotation_version
        create_fake_variants(cls.genome_build)

    def _dump_with_annotated_placeholders(self, output_dir):
        """ Dump external runs, then drop a placeholder annotated VCF next to each sidecar (as the compute
            box + operator would after running VEP and copying results back). """
        annotation_runs = dump_external_annotation_runs(self.variant_annotation_version, output_dir)
        for annotation_run in annotation_runs:
            meta_path = annotation_run.get_dump_metadata_filename(dump_dir=output_dir)
            stem = os.path.basename(meta_path)[:-len(".meta.json")]
            annotated = os.path.join(output_dir, stem + ANNOTATED_VCF_SUFFIX)
            with open(annotated, "w") as f:
                f.write("##fileformat=VCFv4.1\n")
        return annotation_runs

    def test_matching_helpers(self):
        with tempfile.TemporaryDirectory() as output_dir:
            annotation_runs = self._dump_with_annotated_placeholders(output_dir)
            annotation_run = annotation_runs[0]
            meta_path = annotation_run.get_dump_metadata_filename(dump_dir=output_dir)
            meta = parse_dump_metadata(meta_path)

            vav = find_matching_variant_annotation_version(meta)
            self.assertEqual(vav, self.variant_annotation_version)

            matched_run = find_matching_annotation_run(meta, vav)
            self.assertEqual(matched_run, annotation_run)

            self.assertIsNotNone(find_annotated_vcf(meta_path))

    def test_dry_run_matches_without_importing(self):
        with tempfile.TemporaryDirectory() as output_dir:
            annotation_runs = self._dump_with_annotated_placeholders(output_dir)
            report = import_external_annotation_runs(self.genome_build, output_dir, dry_run=True)

            self.assertEqual(len(report["matched"]), len(annotation_runs))
            self.assertFalse(report["imported"])
            self.assertFalse(report["id_mismatch"])
            # Runs untouched - still awaiting external annotation
            for annotation_run in annotation_runs:
                annotation_run.refresh_from_db()
                self.assertEqual(annotation_run.status, AnnotationStatus.EXTERNAL_DUMP_COMPLETED)

    def test_missing_annotated_reported(self):
        with tempfile.TemporaryDirectory() as output_dir:
            # Dump but write NO annotated files - operator hasn't run VEP yet
            dump_external_annotation_runs(self.variant_annotation_version, output_dir)
            report = import_external_annotation_runs(self.genome_build, output_dir, dry_run=True)
            self.assertTrue(report["missing_annotated"])
            self.assertFalse(report["matched"])

    def test_id_mismatch_marks_run_and_continues(self):
        # A genuine divergence is caught by coordinate matching first (matching + §6a both check the range
        # endpoints), so exercise the orchestration's §6a-failure branch by forcing the check to raise: the
        # matched run must be marked ERROR + reported, without aborting the whole import.
        with tempfile.TemporaryDirectory() as output_dir:
            annotation_runs = self._dump_with_annotated_placeholders(output_dir)
            with patch("annotation.external_annotation.verify_annotated_vcf_variant_ids",
                       side_effect=VariantIdAlignmentError("forced divergence")):
                report = import_external_annotation_runs(self.genome_build, output_dir, dry_run=True)

            self.assertEqual(len(report["id_mismatch"]), len(annotation_runs))
            self.assertFalse(report["matched"])
            for annotation_run in annotation_runs:
                annotation_run.refresh_from_db()
                self.assertEqual(annotation_run.status, AnnotationStatus.ERROR)

    def test_import_copies_and_kicks_upload_only(self):
        with tempfile.TemporaryDirectory() as output_dir, tempfile.TemporaryDirectory() as dump_dir:
            annotation_runs = self._dump_with_annotated_placeholders(output_dir)
            with override_settings(ANNOTATION_VCF_DUMP_DIR=dump_dir):
                with patch("annotation.external_annotation.annotation_run_retry") as mock_retry:
                    report = import_external_annotation_runs(self.genome_build, output_dir)

            self.assertEqual(len(report["imported"]), len(annotation_runs))
            self.assertEqual(mock_retry.call_count, len(annotation_runs))
            for _, kwargs in mock_retry.call_args_list:
                self.assertTrue(kwargs["upload_only"])

            for annotation_run in annotation_runs:
                annotation_run.refresh_from_db()
                self.assertIsNotNone(annotation_run.vcf_annotated_filename)
                # Copied into ANNOTATION_VCF_DUMP_DIR
                self.assertEqual(os.path.dirname(annotation_run.vcf_annotated_filename), dump_dir)
                self.assertTrue(os.path.exists(annotation_run.vcf_annotated_filename))
