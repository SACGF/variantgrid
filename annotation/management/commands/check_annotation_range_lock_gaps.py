"""
Report variants left between AnnotationRangeLocks of each build's latest VariantAnnotationVersion, which
the scheduler (searching only above the highest lock) never annotates (#1044). Each gap's pk range is
probed row by row, so a gap holding only other builds' variants is slow - each probe has a time budget and
the gaps that exceed it are listed as unchecked (#2416). Run it on demand, not on deploy.
"""
from django.core.management.base import BaseCommand

from annotation.annotation_versions import get_range_lock_gaps_with_variants
from annotation.models.models_version import VariantAnnotationVersion
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    category = "ops"

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', help="Only check this build (default: every build with annotation)")
        parser.add_argument('--probe-timeout', type=int, default=60,
                            help="Seconds to spend probing each gap before reporting it unchecked (0 = no limit)")

    def handle(self, *args, **options):
        if build_name := options["genome_build"]:
            genome_builds = [GenomeBuild.get_name_or_alias(build_name)]
        else:
            genome_builds = list(GenomeBuild.builds_with_annotation())

        for genome_build in genome_builds:
            vav = VariantAnnotationVersion.latest(genome_build)
            if vav is None:
                self.stdout.write(f"{genome_build}: no VariantAnnotationVersion")
                continue
            gaps = get_range_lock_gaps_with_variants(vav, probe_timeout_seconds=options["probe_timeout"])
            if gaps.with_variants:
                self.stdout.write(self.style.WARNING(
                    f"{vav} has {len(gaps.with_variants)} gap(s) between AnnotationRangeLocks holding variants that "
                    "will never be annotated. Create an AnnotationRangeLock over each gap's variants and the "
                    "annotation scheduler will create its runs. Variant pk ranges:"))
                self._write_ranges(gaps.with_variants)
            if gaps.unchecked:
                self.stdout.write(self.style.WARNING(
                    f"{vav}: {len(gaps.unchecked)} gap(s) not checked, probe exceeded {options['probe_timeout']}s "
                    "(re-run with a larger --probe-timeout). Variant pk ranges:"))
                self._write_ranges(gaps.unchecked)
            if not (gaps.with_variants or gaps.unchecked):
                self.stdout.write(f"{vav}: no gaps between AnnotationRangeLocks hold variants")

    def _write_ranges(self, ranges: list[tuple[int, int]]):
        for min_id, max_id in ranges:
            self.stdout.write(f"  {min_id}-{max_id}")
