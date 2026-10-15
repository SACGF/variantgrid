"""
Backfill ResolvedVariantInfo.p_hgvs / p_hgvs_source / p_hgvs_error (#1383). Only those three columns are written,
so it is cheap and safe to rerun - eg after a new VEP version lands or transcript sequences have been cached.

The default pass never fetches transcript sequences; --retrieve-transcripts is a second, slower pass that fetches
the missing ones from NCBI / Ensembl for rows the first pass couldn't translate (use with --only-missing).
"""
from django.core.management.base import BaseCommand
from django.db.models import Count

from classification.models.classification_variant_info_models import (
    ResolvedVariantInfo,
    p_hgvs_converter,
)
from library.django_utils import chunked_queryset
from snpdb.models.models_genome import GenomeBuild

P_HGVS_FIELDS = ["p_hgvs", "p_hgvs_source", "p_hgvs_error"]


class Command(BaseCommand):
    category = "maintenance"

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', help="Only this build (default all)")
        parser.add_argument('--only-missing', action='store_true', help="Only rows with no p_hgvs")
        parser.add_argument('--retrieve-transcripts', action='store_true',
                            help="Fetch missing transcript sequences (slow, network)")
        parser.add_argument('--chunk-size', type=int, default=1000)

    def handle(self, *args, **options):
        rvi_qs = ResolvedVariantInfo.objects.filter(transcript_version__isnull=False)
        if build_name := options["genome_build"]:
            rvi_qs = rvi_qs.filter(genome_build=GenomeBuild.get_name_or_alias(build_name))
        summary_qs = rvi_qs
        if options["only_missing"]:
            rvi_qs = rvi_qs.filter(p_hgvs__isnull=True)

        num_processed = recalc_p_hgvs_for_queryset(rvi_qs, retrieve_transcripts=options["retrieve_transcripts"],
                                                   chunk_size=options["chunk_size"])
        self.stdout.write(f"Processed {num_processed} ResolvedVariantInfo rows")
        self._print_summary(summary_qs)

    def _print_summary(self, rvi_qs):
        source_labels = dict(ResolvedVariantInfo._meta.get_field("p_hgvs_source").choices)
        self.stdout.write("Rows per p.HGVS source:")
        for row in rvi_qs.values("p_hgvs_source").annotate(count=Count("pk")).order_by("-count"):
            self.stdout.write(f"  {source_labels.get(row['p_hgvs_source'], 'None')}: {row['count']}")

        error_qs = rvi_qs.filter(p_hgvs__isnull=True, p_hgvs_error__isnull=False)
        self.stdout.write(f"Rows with an error: {error_qs.count()}")
        for row in error_qs.values("p_hgvs_error").annotate(count=Count("pk")).order_by("-count")[:10]:
            self.stdout.write(f"  {row['count']}: {row['p_hgvs_error']}")


def recalc_p_hgvs_for_queryset(rvi_qs, retrieve_transcripts: bool = False, chunk_size: int = 1000) -> int:
    rvi_qs = rvi_qs.select_related("genome_build", "variant", "transcript_version",
                                   "allele_info__allele__clingen_allele")
    num_processed = 0
    for chunk_qs in chunked_queryset(rvi_qs, chunk_size):
        for rvi in chunk_qs:
            rvi.set_p_hgvs(converter=p_hgvs_converter(rvi.genome_build, retrieve_transcripts=retrieve_transcripts))
            rvi.save(update_fields=P_HGVS_FIELDS)
            num_processed += 1
    return num_processed
