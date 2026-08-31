import logging
from collections import Counter

from django.core.management.base import BaseCommand
from django.db.models import Q

from annotation.models import (
    AbstractVariantAnnotation,
    VariantAnnotation,
    VariantTranscriptAnnotation,
)
from genes.hgvs import HGVSMatcher
from library.genomics.vcf_enums import VCFSymbolicAllele
from snpdb.models.models_genome import GenomeBuild

PLACEHOLDER_MESSAGES = [
    AbstractVariantAnnotation.SV_HGVS_TOO_LONG_MESSAGE,
    AbstractVariantAnnotation.SV_HGVS_ERROR_MESSAGE,
]
RANGED_SYMBOLIC_ALTS = [VCFSymbolicAllele.DEL, VCFSymbolicAllele.DUP, VCFSymbolicAllele.INV]


class Command(BaseCommand):
    """ Symbolic DEL/DUP/INV are HGVS from coordinates alone (#1571) - rows imported before that
        hold a placeholder message where the string should be, so recalculate them """

    BATCH_SIZE = 1000

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', help="Only process this build (default: all with annotation)")

    def handle(self, *args, **options):
        if genome_build_name := options["genome_build"]:
            genome_builds = [GenomeBuild.get_name_or_alias(genome_build_name)]
        else:
            genome_builds = GenomeBuild.builds_with_annotation()

        for genome_build in genome_builds:
            # We have the transcripts locally - ClinGen would reject these on size anyway
            hgvs_matcher = HGVSMatcher(genome_build, clingen_resolution=False,
                                       allow_alternative_transcript_version=False)
            for klass in [VariantAnnotation, VariantTranscriptAnnotation]:
                self._recalculate(genome_build, hgvs_matcher, klass)

    def _recalculate(self, genome_build: GenomeBuild, hgvs_matcher: HGVSMatcher, klass):
        has_hgvs_g = klass is VariantAnnotation
        fields = ["hgvs_g", "hgvs_c"] if has_hgvs_g else ["hgvs_c"]

        placeholder_q = Q(hgvs_c__in=PLACEHOLDER_MESSAGES)
        if has_hgvs_g:
            placeholder_q |= Q(hgvs_g__in=PLACEHOLDER_MESSAGES)

        qs = klass.objects.filter(placeholder_q,
                                  version__genome_build=genome_build,
                                  variant__alt__seq__in=RANGED_SYMBOLIC_ALTS)
        qs = qs.select_related("variant__locus__contig", "variant__locus__ref", "variant__alt",
                               "transcript_version__transcript")

        results = Counter()
        records = []
        for record in qs.iterator():
            variant_coordinate = record.variant.coordinate
            if has_hgvs_g and record.hgvs_g in PLACEHOLDER_MESSAGES:
                try:
                    record.hgvs_g = hgvs_matcher.variant_coordinate_to_g_hgvs(variant_coordinate)
                    results["hgvs_g ok"] += 1
                except Exception as e:
                    logging.warning("Error calculating g.HGVS for '%s': %s", variant_coordinate.format_short(), e)
                    record.hgvs_g = VariantAnnotation.SV_HGVS_ERROR_MESSAGE
                    results["hgvs_g error"] += 1

            if record.hgvs_c in PLACEHOLDER_MESSAGES:
                hgvs_c = None  # c.HGVS is ok to be blank - the g.HGVS carries the answer
                if transcript_version := record.transcript_version:
                    try:
                        hgvs_c = str(hgvs_matcher.variant_coordinate_to_hgvs_variant(variant_coordinate,
                                                                                     transcript_version.accession))
                        results["hgvs_c ok"] += 1
                    except Exception as e:
                        logging.debug("Error calculating c.HGVS for '%s'/'%s': %s",
                                      variant_coordinate.format_short(), transcript_version.accession, e)
                        results["hgvs_c error"] += 1
                record.hgvs_c = hgvs_c

            records.append(record)
            if len(records) >= self.BATCH_SIZE:
                klass.objects.bulk_update(records, fields, batch_size=self.BATCH_SIZE)
                records = []

        if records:
            klass.objects.bulk_update(records, fields, batch_size=self.BATCH_SIZE)

        logging.info("%s/%s: %s", genome_build, klass.__name__, dict(results))
