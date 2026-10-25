"""
Switch on UTRAnnotator / AVI / SpliceVault (#2426) and AlphaGenome splicing (#2433) for the VariantAnnotationVersions
in use, after setting their files in settings.ANNOTATION[build]["vep_config"]: stamps the pins, flags versions for the
backfill and prints its steps. @see annotation/vep_plugin_backfill.py
"""
from django.core.management.base import BaseCommand, CommandError

from annotation.enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vep_annotation import get_missing_vep_files
from annotation.vep_config import disallowed_non_commercial_keys
from annotation.vep_plugin_backfill import (
    ACTIVE_STATUSES,
    SOURCE_PINS,
    backfill_steps,
    stamp_plugin_pins,
)
from snpdb.models.models_genome import GenomeBuild

# Runs whose VEP ran (or is running) on the old command - importing them once the inserter expects the new CSQ
# fields fails with VEPMissingColumnsError
_VEP_RUN_NOT_IMPORTED = [AnnotationStatus.EXTERNAL_DUMP_COMPLETED, AnnotationStatus.ANNOTATION_STARTED,
                         AnnotationStatus.ANNOTATION_COMPLETED, AnnotationStatus.UPLOAD_STARTED]


class Command(BaseCommand):
    category = "ops"
    help = ("Stamp UTRAnnotator / AVI / SpliceVault / AlphaGenome splicing pins from settings on the versions in use "
            "and print the backfill")

    def add_arguments(self, parser):
        parser.add_argument("--genome-build", help="Default: every build with annotation")
        parser.add_argument("--dry-run", action="store_true", help="Say what would be stamped, changing nothing")
        parser.add_argument("--force", action="store_true",
                            help="Stamp even with runs VEP-annotated but not yet imported (they will fail to "
                                 "import and need re-running)")

    def handle(self, *args, **options):
        if build_name := options["genome_build"]:
            genome_builds = [GenomeBuild.get_name_or_alias(build_name)]
        else:
            genome_builds = list(GenomeBuild.builds_with_annotation())

        problems = []
        for genome_build in genome_builds:
            for key, filename in get_missing_vep_files(genome_build).items():
                problems.append(f"{genome_build} {key}: {filename} MISSING")
        if problems:
            raise CommandError("VEP would fail on missing files - install them or unset them in settings:\n"
                               + "\n".join(problems))

        for genome_build in genome_builds:
            if keys := disallowed_non_commercial_keys(genome_build.settings["vep_config"]):
                self.stdout.write(f"{genome_build}: not stamping {', '.join(keys)} - non-commercial data, and "
                                  f"ANNOTATION_NON_COMMERCIAL_DATA_ALLOWED is False")

        if not options["force"]:
            pending_runs = AnnotationRun.objects.filter(
                annotation_range_lock__version__genome_build__in=genome_builds,
                annotation_range_lock__version__status__in=ACTIVE_STATUSES,
                pipeline_type=VariantAnnotationPipelineType.STANDARD,
                status__in=_VEP_RUN_NOT_IMPORTED)
            if count := pending_runs.count():
                raise CommandError(f"{count} STANDARD annotation runs have VEP output from the current command that "
                                   f"isn't imported yet, and would fail to import without the new CSQ fields. Let "
                                   f"the import lane drain (or re-run them afterwards) and use --force")

        stamps = []
        for genome_build in genome_builds:
            stamps.extend(stamp_plugin_pins(VariantAnnotationVersion, AnnotationRangeLock,
                                            genome_build_name=genome_build.name, dry_run=options["dry_run"]))

        verb = "Would stamp" if options["dry_run"] else "Stamped"
        for stamp in stamps:
            vav = stamp.variant_annotation_version
            if stamp.changed_pins:
                changes = ", ".join(f"{pin}: {old} -> {new}" for pin, (old, new) in stamp.changed_pins.items())
                self.stdout.write(f"{verb} {vav}: {changes}")
            else:
                pins = ", ".join(f"{pin}={getattr(vav, pin)}" for pin in SOURCE_PINS)
                self.stdout.write(f"{vav}: already matches settings ({pins})")

        if any(stamp.changed_pins for stamp in stamps) and not options["dry_run"]:
            self.stdout.write("\nStart annotation_workers (stop them before running this, as get_vep_version is "
                              "cached per process) and restart the other celery workers and gunicorn, which read "
                              "settings at start. A VEP run on the old settings latches VEPVersionMismatchError on "
                              "its version - clear it with the re-check on the annotation runs page.")

        if to_backfill := [stamp.variant_annotation_version for stamp in stamps if stamp.needs_backfill]:
            self.stdout.write("\nThen backfill the variants annotated before they were switched on, running each "
                              "build's commands in order (https://github.com/SACGF/variantgrid/issues/2436):")
            for vav in to_backfill:
                self.stdout.write("")
                self.stdout.write("\n".join(backfill_steps(vav, AnnotationRangeLock)))
