#!/usr/bin/env python3

from django.core.management.base import BaseCommand
import logging

from genes.canonical_transcripts.create_canonical_transcripts import create_canonical_transcript_collection
from library.utils import invert_dict
from snpdb.models import GenomeBuild, AnnotationConsortium


class Command(BaseCommand):

    def add_arguments(self, parser):
        consortia = [ac[1] for ac in AnnotationConsortium.choices]
        builds = [gb.name for gb in GenomeBuild.builds_with_annotation()]

        parser.add_argument('--genome-build', choices=builds, required=True)
        parser.add_argument('--annotation-consortium', choices=consortia, required=True)
        parser.add_argument('filename')

    def handle(self, *args, **options):
        annotation_consortium_name = options["annotation_consortium"]
        build_name = options["genome_build"]
        filename = options['filename']

        ac_dict = invert_dict(dict(AnnotationConsortium.choices))
        annotation_consortium = ac_dict[annotation_consortium_name]
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        collection = create_canonical_transcript_collection(genome_build, annotation_consortium, filename)

        logging.info("Finished - created as %s - please link EnrichmentKit to this via admin tool", collection.pk)
