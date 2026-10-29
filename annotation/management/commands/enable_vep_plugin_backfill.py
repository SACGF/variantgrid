"""
Switch on a VEP plugin / custom VCF (those in vep_plugin_backfill.SOURCE_PINS, first #2426 and #2433) for the
VariantAnnotationVersions in use, after setting its file in settings.ANNOTATION[build]["vep_config"]: stamps the pins,
flags versions for the backfill and prints its steps. @see annotation/vep_plugin_backfill.py
"""
import os
import stat

from django.core.management.base import BaseCommand, CommandError

from annotation.enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.models.models_annotation_run import AnnotationRangeLock, AnnotationRun
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vep_annotation import get_missing_vep_files
from annotation.vep_config import disallowed_non_commercial_keys
from annotation.vep_plugin_backfill import (
    ACTIVE_STATUSES,
    SOURCE_PINS,
    backfill_script,
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
    help = ("Switch on VEP plugin / custom VCF files newly set in settings for the VariantAnnotationVersions in use, "
            "without a new version: stamps their pins and prints the backfill for variants already annotated, "
            "or writes it as a bash script with --bash-script")

    def add_arguments(self, parser):
        parser.add_argument("--genome-build", help="Default: every build with annotation")
        parser.add_argument("--dry-run", action="store_true", help="Say what would be stamped, changing nothing")
        parser.add_argument("--force", action="store_true",
                            help="Stamp even with runs VEP-annotated but not yet imported (they will fail to "
                                 "import and need re-running)")
        parser.add_argument("--bash-script", metavar="PATH",
                            help="Also write the backfill steps as a bash script to PATH, to run under screen "
                                 "or nohup (with --dry-run, the script for what would be flagged)")

    @staticmethod
    def _stamp(genome_builds, dry_run: bool):
        stamps = []
        for genome_build in genome_builds:
            stamps.extend(stamp_plugin_pins(VariantAnnotationVersion, AnnotationRangeLock,
                                            genome_build_name=genome_build.name, dry_run=dry_run))
        return stamps

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

        stamps = self._stamp(genome_builds, dry_run=True)
        # Only a pin change puts runs at risk: a re-run for a later backfill pass stamps nothing and can go ahead
        # with the annotation lane busy
        if any(stamp.changed_pins for stamp in stamps) and not options["force"]:
            pending_runs = AnnotationRun.objects.filter(
                annotation_range_lock__version__genome_build__in=genome_builds,
                annotation_range_lock__version__status__in=ACTIVE_STATUSES,
                pipeline_type=VariantAnnotationPipelineType.STANDARD,
                status__in=_VEP_RUN_NOT_IMPORTED)
            if count := pending_runs.count():
                raise CommandError(f"{count} STANDARD annotation runs have VEP output from the current command that "
                                   f"isn't imported yet, and would fail to import without the new CSQ fields. Let "
                                   f"the import lane drain (or re-run them afterwards) and use --force")
        if not options["dry_run"]:
            stamps = self._stamp(genome_builds, dry_run=False)

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

        to_backfill = [stamp.variant_annotation_version for stamp in stamps if stamp.needs_backfill]
        if to_backfill:
            self.stdout.write("\nThen backfill the variants annotated before they were switched on, running each "
                              "build's commands in order "
                              "(https://github.com/SACGF/variantgrid/issues/2433#issuecomment-6393595562), or "
                              "re-run this with --bash-script PATH to write them as one script:")
            for vav in to_backfill:
                self.stdout.write("")
                self.stdout.write("\n".join(backfill_steps(vav, AnnotationRangeLock)))

        if script_path := options["bash_script"]:
            with open(script_path, "w", encoding="utf-8") as f:
                f.write(backfill_script(to_backfill, AnnotationRangeLock))
            os.chmod(script_path, os.stat(script_path).st_mode | stat.S_IXUSR)
            self.stdout.write(f"\nWrote the backfill as a script to {script_path}")
