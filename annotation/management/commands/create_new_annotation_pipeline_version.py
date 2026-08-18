#!/usr/bin/env python3
"""
Register the installed version of a non-VEP annotation pipeline's tool (#720).

The counterpart of create_new_variant_annotation_version, and the way an upgrade is rolled out: install
the new binary/bundle, run this, then promote the NEW AnnotationPipelineVersion on the annotation runs
page. Promoting is what makes the scheduler create runs against it, so every range is re-annotated by the
machinery that already backfills a newly-enabled pipeline.

The first version for a build is registered ACTIVE - there is no earlier one for it to supersede.
"""
import logging

from django.core.management.base import BaseCommand, CommandError

from annotation.models.models_enums import VariantAnnotationPipelineType
from annotation.pipelines import PIPELINES, get_runner, versioned_pipeline_types
from snpdb.models.models_genome import GenomeBuild

class Command(BaseCommand):
    def add_arguments(self, parser):
        parser.add_argument("--pipeline-type", required=True,
                            choices=[pt.value for pt in versioned_pipeline_types()],
                            help="AnnotationRun pipeline_type, eg 'A' for AnnotSV")
        parser.add_argument("--genome-build")

    def handle(self, *args, **options):
        pipeline_type = VariantAnnotationPipelineType(options["pipeline_type"])
        if not PIPELINES[pipeline_type].enabled:
            raise CommandError(f"{pipeline_type.label} is not enabled in settings")
        if build_name := options.get("genome_build"):
            genome_builds = [GenomeBuild.get_name_or_alias(build_name)]
        else:
            genome_builds = list(GenomeBuild.builds_with_annotation())

        runner = get_runner(pipeline_type)
        for genome_build in genome_builds:
            pipeline_version, created = runner.get_or_create_current_version(genome_build)
            if created:
                logging.info("Created: %s", pipeline_version)
                if pipeline_version.is_active:
                    self.stdout.write(self.style.SUCCESS(f"{genome_build}: registered {pipeline_version} "
                                                         f"as ACTIVE - the scheduler will annotate against it"))
                else:
                    self.stdout.write(self.style.SUCCESS(
                        f"{genome_build}: registered {pipeline_version} as NEW. Promote it on the "
                        f"annotation runs page to re-annotate every range with it."))
            else:
                self.stdout.write(f"{genome_build}: already registered - {pipeline_version} "
                                  f"[{pipeline_version.get_status_display()}]")
