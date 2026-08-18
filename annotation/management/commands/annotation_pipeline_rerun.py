#!/usr/bin/env python3
"""
Re-run a non-VEP annotation pipeline over ranges it has already finished.

Two uses, both of which the pipeline's own AnnotationRuns make possible (#720):

* the tool's version rolled - `--stale` picks the FINISHED runs that used an older
  AnnotationPipelineVersion than the one installed now;
* the pipeline was enabled after those ranges were annotated and the runs finished empty.

Safe to run on FINISHED runs: a supplementary run owns no VariantAnnotation rows (the annotation_run FK
points at the VEP run that created them), so resetting one throws away nothing but the run's own state.
"""
import logging

from django.core.management.base import BaseCommand, CommandError

from annotation.models.models import AnnotationRun
from annotation.models.models_enums import AnnotationStatus, VariantAnnotationPipelineType
from annotation.pipelines import PIPELINES
from annotation.pipelines.annotsv import get_current_annotsv_version
from annotation.tasks.annotate_variants import reset_annotation_run_for_retry
from annotation.tasks.annotation_scheduler_task import dispatch_annotation_runs
from snpdb.models import GenomeBuild

# How to resolve "the version the deployment would use right now", per pipeline type.
CURRENT_VERSION_RESOLVERS = {
    VariantAnnotationPipelineType.ANNOTSV: get_current_annotsv_version,
}


class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("--pipeline-type", required=True,
                            choices=[pt.value for pt in CURRENT_VERSION_RESOLVERS],
                            help="AnnotationRun pipeline_type, eg 'A' for AnnotSV")
        parser.add_argument("--genome-build", required=True)
        parser.add_argument("--stale", action="store_true",
                            help="Only runs whose pipeline_version isn't the one installed now")
        parser.add_argument("--limit", type=int, help="Reset at most this many runs")
        parser.add_argument("--dry-run", action="store_true")

    def handle(self, *args, **options):
        pipeline_type = VariantAnnotationPipelineType(options["pipeline_type"])
        if not PIPELINES[pipeline_type].enabled:
            raise CommandError(f"{pipeline_type.label} is not enabled in settings")
        genome_build = GenomeBuild.get_name_or_alias(options["genome_build"])

        qs = AnnotationRun.objects.filter(
            pipeline_type=pipeline_type,
            status=AnnotationStatus.FINISHED,
            annotation_range_lock__version__genome_build=genome_build,
        )
        if options["stale"]:
            current = CURRENT_VERSION_RESOLVERS[pipeline_type](genome_build)
            qs = qs.exclude(pipeline_version=current)
            self.stdout.write(f"Current version: {current}")
        qs = qs.order_by("annotation_range_lock__min_variant_id")
        if limit := options["limit"]:
            qs = qs[:limit]

        run_ids = list(qs.values_list("pk", flat=True))
        self.stdout.write(f"{len(run_ids)} {pipeline_type.label} run(s) to re-run")
        if options["dry_run"] or not run_ids:
            return

        for run_id in run_ids:
            reset_annotation_run_for_retry.apply_async((run_id,))
        logging.info("Queued %d resets", len(run_ids))
        # One kick is enough - each reset triggers the dispatcher again as it completes.
        version_id = AnnotationRun.objects.get(pk=run_ids[0]).annotation_range_lock.version_id
        dispatch_annotation_runs.si(version_id).apply_async()
        self.stdout.write(self.style.SUCCESS(f"Queued {len(run_ids)} run(s) for re-annotation"))
