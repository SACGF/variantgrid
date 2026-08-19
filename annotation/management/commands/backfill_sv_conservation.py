#!/usr/bin/env python3
"""
Recompute the structural-variant conservation (phastCons/phyloP _max) columns with pyBigWig (#1657).

Two things put wrong values into those columns on systems that annotated SVs with the pyBigWig stage:

 - the scoring window started one base early, including the symbolic ALT's anchor base, which VEP
   excludes. Only bites when that base holds the span max, but then it can be metres out (a phyloP
   7.75 where VEP says 3.93).
 - a run whose sidecar stage failed imported nulls, which read as "no conservation here".

Both are corrected in place by re-scoring each SV against the build's bigWig tracks: cheap (ms/SV),
and far lighter than re-running the annotation. Only runs flagged sv_conservation_pybigwig are touched -
a run VEP scored itself already holds the values VEP would give - and only where the stored value really
differs, so re-running is a no-op. Versions re-scored end up with backfilled_sv_conservation set.
"""
import logging
from collections import Counter

from django.core.management.base import BaseCommand

from annotation.models import (
    AnnotationRun,
    AnnotationStatus,
    VariantAnnotation,
    VariantAnnotationPipelineType,
    VariantAnnotationVersion,
)
from annotation.sv_conservation import (
    ConservationTrack,
    get_sv_conservation_tracks,
    score_sv_variants,
)
from snpdb.models.models_genome import GenomeBuild

# VEP recorded these as 3-decimal text, so a re-score of a VEP-computed value lands ~1e-7 away from
# what's stored. Only a difference well above that is a value worth rewriting - the smallest real
# off-by-one seen in production data moved a value by 0.003.
MIN_CHANGE = 1e-4


class Command(BaseCommand):
    help = "Recompute SV conservation (phastCons/phyloP max) columns with pyBigWig (#1657)"

    def add_arguments(self, parser):
        parser.add_argument("--genome-build", help="Default: every build with annotation")
        parser.add_argument("--batch-size", type=int, default=1000,
                            help="Variants scored (and rows updated) per batch")
        parser.add_argument("--threads", type=int, default=None,
                            help="Thread pool size (default: settings.ANNOTATION_VEP_FORK)")
        parser.add_argument("--dry-run", action="store_true", help="Report what would change, write nothing")

    def handle(self, *args, **options):
        if genome_build_name := options["genome_build"]:
            genome_builds = [GenomeBuild.get_name_or_alias(genome_build_name)]
        else:
            genome_builds = list(GenomeBuild.builds_with_annotation())

        counts = Counter()
        for genome_build in genome_builds:
            self._handle_genome_build(genome_build, counts, options)

        verb = "would update" if options["dry_run"] else "updated"
        self.stdout.write(f"Scanned {counts['runs']} pyBigWig-scored SV annotation runs / {counts['rows']} rows: "
                          f"{verb} {counts['rows_changed']} rows "
                          f"({counts['values_corrected']} corrected, {counts['values_filled']} filled in)")
        self.stdout.write(f"Marked {counts['versions_backfilled']} annotation version(s) backfilled; "
                          f"skipped {counts['runs_no_tracks']} runs on builds with no conservation "
                          f"tracks configured")

    def _handle_genome_build(self, genome_build: GenomeBuild, counts: Counter, options: dict):
        run_qs = AnnotationRun.objects.filter(
            pipeline_type=VariantAnnotationPipelineType.STRUCTURAL_VARIANT,
            status=AnnotationStatus.FINISHED,
            annotated_count__gt=0,
            annotation_range_lock__version__genome_build=genome_build,
            sv_conservation_pybigwig=True,
        ).select_related("annotation_range_lock__version").order_by("pk")

        tracks_by_version_id: dict[int, list[ConservationTrack]] = {}
        rescored_version_ids = set()
        total = run_qs.count()
        logging.info("%s: %d finished pyBigWig-scored SV annotation runs", genome_build, total)
        for i, annotation_run in enumerate(run_qs, start=1):
            version = annotation_run.annotation_range_lock.version
            if version.pk not in tracks_by_version_id:
                # Which conservation columns this annotation version has, rather than what the current
                # install would annotate - an older version may predate some of the tracks.
                tracks_by_version_id[version.pk] = get_sv_conservation_tracks(
                    genome_build, columns_version=version.columns_version, vep_version=version.vep)
            if tracks := tracks_by_version_id[version.pk]:
                counts["runs"] += 1
                rescored_version_ids.add(version.pk)
                self._handle_run(annotation_run, version, tracks, counts, options)
            else:
                counts["runs_no_tracks"] += 1
            if i % 100 == 0:
                logging.info("%s: processed %d of %d runs", genome_build, i, total)

        # Only versions every run of which we could actually re-score are clean - a build whose bigWig
        # data files aren't installed has to stay flagged so the manual task comes back
        if rescored_version_ids and not options["dry_run"]:
            counts["versions_backfilled"] += VariantAnnotationVersion.objects.filter(
                pk__in=rescored_version_ids).update(backfilled_sv_conservation=True)

    def _handle_run(self, annotation_run: AnnotationRun, version, tracks: list[ConservationTrack],
                    counts: Counter, options: dict):
        columns = [t.db_column for t in tracks]
        # version as well as the run so this prunes to the one partition rather than hitting every
        # version's annotation_run index
        qs = VariantAnnotation.objects.filter(version=version, annotation_run=annotation_run) \
            .select_related("variant__locus__contig").only("pk", "variant", *columns)

        batch: list[VariantAnnotation] = []
        for variant_annotation in qs.iterator(chunk_size=options["batch_size"]):
            batch.append(variant_annotation)
            if len(batch) >= options["batch_size"]:
                self._handle_batch(batch, tracks, columns, counts, options)
                batch = []
        if batch:
            self._handle_batch(batch, tracks, columns, counts, options)

    def _handle_batch(self, batch: list[VariantAnnotation], tracks: list[ConservationTrack],
                      columns: list[str], counts: Counter, options: dict):
        counts["rows"] += len(batch)
        variants = []
        for variant_annotation in batch:
            variant = variant_annotation.variant
            variants.append((variant.pk, variant.locus.contig.name, variant.locus.position,
                             variant.end, variant.svlen))
        scored = score_sv_variants(variants, tracks, threads=options["threads"])

        to_update = []
        for variant_annotation in batch:
            values = scored.get(variant_annotation.variant_id, {})
            changed = False
            for column in columns:
                if (value := values.get(column)) is None:
                    continue  # bigWig has no data over the span - leave whatever is there
                existing = getattr(variant_annotation, column)
                if existing is None:
                    counts["values_filled"] += 1
                elif abs(existing - value) > MIN_CHANGE:
                    counts["values_corrected"] += 1
                else:
                    continue
                setattr(variant_annotation, column, value)
                changed = True
            if changed:
                to_update.append(variant_annotation)

        counts["rows_changed"] += len(to_update)
        if to_update and not options["dry_run"]:
            VariantAnnotation.objects.bulk_update(to_update, columns, batch_size=options["batch_size"])
