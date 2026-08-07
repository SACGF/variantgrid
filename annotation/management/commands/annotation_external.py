#!/usr/bin/env python3
"""
External annotation runs (#1568): single entry point for external annotation.

The same command both dumps (write VCFs + metadata, park runs awaiting external annotation) and imports
(re-import annotated VCFs), selected by --dump / --import. Run against a NEW (not yet ACTIVE)
VariantAnnotationVersion: the normal scheduler only operates on the latest ACTIVE version, so it will not
touch a NEW version's range locks.

See claude/plans/1568_external_annotation_runs_plan.md.
"""
import os

from django.core.management.base import BaseCommand, CommandError

from annotation.external_annotation import dump_external_annotation_runs, dump_existing_annotation_runs, \
    import_external_annotation_runs, _require_sv_offload_supported, DEFAULT_MIN_EXTERNAL_VARIANTS
from annotation.models.models import VariantAnnotationVersion
from annotation.models.models_enums import VariantAnnotationPipelineType
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    help = "Dump/import external annotation runs (#1568)"

    def add_arguments(self, parser):
        mode = parser.add_mutually_exclusive_group(required=True)
        mode.add_argument("--dump", action="store_true",
                          help="Create + dump all runs for the NEW version, parking them awaiting external annotation")
        mode.add_argument("--dump-existing", action="store_true", dest="dump_existing",
                          help="Adopt existing CREATED runs (already scheduled but not yet annotated), marking "
                               "them external and dumping them; use --leave to keep some on the local pipeline")
        mode.add_argument("--import", action="store_true", dest="import_mode",
                          help="Re-import annotated VCFs, matching them back to local runs")

        parser.add_argument("--genome-build", required=True)
        parser.add_argument("--pipeline-type", default=VariantAnnotationPipelineType.STANDARD,
                            choices=[VariantAnnotationPipelineType.STANDARD,
                                     VariantAnnotationPipelineType.STRUCTURAL_VARIANT],
                            help="STANDARD (small variant), or STRUCTURAL_VARIANT (SV, VEP-only offload; "
                                 "requires AnnotSV disabled). --dump-existing is the recommended SV path - it "
                                 "adopts the scheduler's SV-pipeline runs rather than dumping mostly-empty locks")
        parser.add_argument("--vav-status", default=VariantAnnotationVersion.Status.NEW,
                            choices=[VariantAnnotationVersion.Status.NEW, VariantAnnotationVersion.Status.ACTIVE],
                            help="--dump-existing: status of the VariantAnnotationVersion whose CREATED runs to "
                                 "adopt (default NEW - runs created by 'Run scheduler against NEW', which the "
                                 "dispatcher won't touch; use ACTIVE to offload the version being annotated)")
        parser.add_argument("--leave", type=int, default=0,
                            help="--dump-existing: leave this many lowest-id CREATED runs on the local pipeline "
                                 "(so annotation can be parallelised across machines)")
        parser.add_argument("--min-variants", type=int, default=DEFAULT_MIN_EXTERNAL_VARIANTS,
                            help="--dump/--dump-existing: keep runs with fewer than this many variants on the "
                                 "local pipeline instead of parking them external (default %(default)s)")
        parser.add_argument("--output-dir", help="--dump/--dump-existing: directory to write VCFs + metadata")
        parser.add_argument("--input-dir", help="--import: directory of annotated VCFs + sidecar metadata")
        parser.add_argument("--dry-run", action="store_true",
                            help="--import: list matches/unmatched without importing")

    def handle(self, *args, **options):
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])
        try:
            _require_sv_offload_supported(options["pipeline_type"])
        except ValueError as e:
            raise CommandError(str(e))
        if options["dump"]:
            self._run_dump(genome_build, options)
        elif options["dump_existing"]:
            self._run_dump_existing(genome_build, options)
        else:
            self._run_import(genome_build, options)

    def _run_dump(self, genome_build, options):
        output_dir = options.get("output_dir")
        if not output_dir:
            raise CommandError("--dump requires --output-dir")

        pipeline_type = options["pipeline_type"]
        variant_annotation_version = VariantAnnotationVersion.latest(
            genome_build, status=VariantAnnotationVersion.Status.NEW)
        if variant_annotation_version is None:
            raise CommandError(
                f"No NEW VariantAnnotationVersion for {genome_build} - create one with "
                f"create_new_variant_annotation_version and leave it NEW until external annotation completes")

        annotation_runs = dump_external_annotation_runs(variant_annotation_version, output_dir,
                                                        pipeline_type=pipeline_type,
                                                        min_variants=options["min_variants"])
        self.stdout.write(
            f"Dumped {len(annotation_runs)} external annotation run(s) for {variant_annotation_version} "
            f"into {output_dir} (runs with < {options['min_variants']} variants kept on the local pipeline)")

    def _run_dump_existing(self, genome_build, options):
        output_dir = options.get("output_dir")
        if not output_dir:
            raise CommandError("--dump-existing requires --output-dir")

        leave = options["leave"]
        if leave < 0:
            raise CommandError("--leave must be >= 0")

        pipeline_type = options["pipeline_type"]
        vav_status = options["vav_status"]
        variant_annotation_version = VariantAnnotationVersion.latest(genome_build, status=vav_status)
        if variant_annotation_version is None:
            raise CommandError(f"No {vav_status} VariantAnnotationVersion for {genome_build}")

        annotation_runs = dump_existing_annotation_runs(variant_annotation_version, output_dir,
                                                        pipeline_type=pipeline_type, leave=leave,
                                                        min_variants=options["min_variants"])
        self.stdout.write(
            f"Dumped {len(annotation_runs)} existing CREATED run(s) as external for "
            f"{variant_annotation_version} into {output_dir} (left {leave} on the local pipeline; runs with "
            f"< {options['min_variants']} variants kept on the local pipeline)")

    def _run_import(self, genome_build, options):
        input_dir = options.get("input_dir")
        if not input_dir:
            raise CommandError("--import requires --input-dir")
        if not os.path.isdir(input_dir):
            raise CommandError(f"--input-dir {input_dir!r} is not a directory")

        pipeline_type = options["pipeline_type"]
        dry_run = options["dry_run"]

        def emit(category, message):
            # Stream each file's outcome + reason as it's decided, rather than batching to the end of a
            # potentially very long (~1k file) run.
            if category == "matched":
                self.stdout.write(f"[dry-run match] {message}")
            elif category == "imported":
                self.stdout.write(f"[import] {message}")
            elif category in ("missing_annotated", "unmatched"):
                self.stdout.write(self.style.WARNING(f"[skip] {message}"))
            elif category == "id_mismatch":
                self.stdout.write(self.style.ERROR(f"[error] {message}"))
            self.stdout.flush()

        report = import_external_annotation_runs(genome_build, input_dir,
                                                 pipeline_type=pipeline_type, dry_run=dry_run, emit=emit)

        verb = "would import" if dry_run else "imported"
        self.stdout.write(
            f"External annotation import: {verb} {len(report['matched']) + len(report['imported'])}, "
            f"skipped {len(report['unmatched']) + len(report['missing_annotated'])}, "
            f"id-mismatch {len(report['id_mismatch'])}")
