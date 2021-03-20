#!/usr/bin/env python3

from django.core.management.base import BaseCommand
from django.db.models import Q, Func, Value, F
import logging

from annotation.models import VariantAnnotation, VariantTranscriptAnnotation, VariantAnnotationVersion, \
    TranscriptVersion, defaultdict
from snpdb.models.models_genome import GenomeBuild


class Command(BaseCommand):
    """ This should only need to be run on legacy data, with variant annotations that were run before the transcript
        versions it used were inserted. Now we ensure the gene_annotation_release is set so this shouldn't happen """

    def handle(self, *args, **options):
        for genome_build in GenomeBuild.builds_with_annotation():
            for vav in VariantAnnotationVersion.objects.filter(genome_build=genome_build):
                self.fix_variant_annotation_version(vav)

    @staticmethod
    def fix_variant_annotation_version(vav: VariantAnnotationVersion):
        t_qs = TranscriptVersion.objects.filter(transcript__annotation_consortium=vav.annotation_consortium,
                                                genome_build=vav.genome_build)
        transcript_versions_by_id = defaultdict(dict)
        for pk, transcript_id, version in t_qs.values_list("pk", "transcript_id", "version"):
            transcript_versions_by_id[transcript_id][version] = pk

        for klass in [VariantAnnotation, VariantTranscriptAnnotation]:
            missing_qs = klass.objects.filter(Q(transcript__isnull=True) | Q(transcript_version__isnull=True),
                                              version=vav, hgvs_c__isnull=False)
            split_func = Func(F("hgvs_c"), Value(":"), Value(1), function="split_part")
            records = []
            for pk, feature in missing_qs.annotate(feature=split_func).values_list("pk", "feature"):
                t_id, version = TranscriptVersion.get_transcript_id_and_version(feature)
                transcript_versions = transcript_versions_by_id.get(t_id)
                if transcript_versions:
                    transcript_id = t_id  # Know it's valid to link
                    transcript_version_id = transcript_versions.get(version)
                    if transcript_version_id is None:
                        logging.warning(f"Have transcript '{transcript_id}' but no version: '{version}'")
                    records.append(klass(pk=pk, transcript_id=transcript_id, transcript_version_id=transcript_version_id))

            print(f"Updating {len(records)} {klass} records")
            if records:
                klass.objects.bulk_update(records, fields=["transcript_id", "transcript_version_id"], batch_size=2000)
