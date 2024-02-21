from django.core.management import BaseCommand

from annotation.models import VariantAnnotation, VariantTranscriptAnnotation, VariantAnnotationPipelineType, \
    AbstractVariantAnnotation
from genes.hgvs import HGVSMatcher, HGVSException
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        Ensembl VEP doesn't populate HGVS for SVs
        @see https://github.com/Ensembl/ensembl-vep/issues/1222#issuecomment-1188788614
    """
    @staticmethod
    def _calculate_hgvs(matcher, qs):
        """ This can often be slow, eg 1/second """
        va: AbstractVariantAnnotation
        for va in qs:
            try:
                variant_coordinate = va.variant.coordinate
                if va.transcript_version:
                    transcript_accession = va.transcript_version.accession
                    va.hgvs_c = matcher.variant_coordinate_to_hgvs_variant(variant_coordinate, transcript_accession)
                    yield va
            except (ValueError, HGVSException):
                pass

    def _update_annotation(self, genome_build, matcher, klass):
        print(f"Setting hgvs_c on {genome_build.name} / {klass}")

        va_qs = klass.objects.filter(version__genome_build=genome_build,
                                     annotation_run__pipeline_type=VariantAnnotationPipelineType.CNV,
                                     variant__alt__seq__startswith='<')
        va_qs.select_related("variant", "variant__locus", "transcript_version")
        klass.objects.bulk_update(self._calculate_hgvs(matcher, va_qs),
                                  fields=["hgvs_c"], batch_size=50)

    def handle(self, *args, **options):

        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build)
            self._update_annotation(genome_build, matcher, VariantAnnotation)
            self._update_annotation(genome_build, matcher, VariantTranscriptAnnotation)
