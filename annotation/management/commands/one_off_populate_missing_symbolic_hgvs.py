import time

from django.core.management import BaseCommand

from annotation.models import VariantAnnotation, VariantTranscriptAnnotation
from genes.hgvs import HGVSMatcher, HGVSException
from snpdb.models import GenomeBuild, Variant, VariantCoordinate


class Command(BaseCommand):
    """
        Ensembl VEP doesn't populate HGVS for SVs
        @see https://github.com/Ensembl/ensembl-vep/issues/1222#issuecomment-1188788614

        The thing that made this slow was retrieving the ref/alt from genome sequences, so change to do it
        By variant so we only need to do that once

    """
    def _update_annotation(self, v: Variant, variant_coordinate: VariantCoordinate, matcher, records, klass):
        va_list = []
        for va in klass.objects.filter(variant=v, hgvs_c__isnull=True):
            if va.transcript_version:
                # This is not super perfect - as normalization may affect position, but we want a quick reject
                if va.transcript_version.start > variant_coordinate.position:
                    continue
                elif variant_coordinate.end > va.transcript_version.end + 1:
                    continue

                transcript_accession = va.transcript_version.accession
                try:
                    va.hgvs_c = matcher.variant_coordinate_to_hgvs_variant(variant_coordinate, transcript_accession)
                    va_list.append(va)
                except (ValueError, HGVSException) as e:
                    # print(f"FAILED: {transcript_accession}: {e} - {quick_reject=}")
                    pass


        records.extend(va_list)
        self._bulk_update(klass, records)

    def _bulk_update(self, klass, records, min_batch_size=1000):
        if len(records) >= min_batch_size:
            klass.objects.bulk_update(records, fields=["hgvs_c"])
            records.clear()

    def handle(self, *args, **options):
        MAX_SIZE = 100_000

        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build)
            qs = Variant.objects.filter(Variant.get_contigs_q(genome_build))
            symbolic_qs = qs.filter(alt__seq__startswith='<')
            total = symbolic_qs.count()
            print(f"{genome_build.name=} has {total} variants to update annotation")

            va_list = []
            vta_list = []

            last_update = time.time()
            for i, v in enumerate(symbolic_qs):
                # Do this lookup of ref/alt once as it's expensive..
                variant_coordinate = v.coordinate
                if abs(variant_coordinate.svlen) > MAX_SIZE:
                    print(f"Skipping calculating HGVS for {variant_coordinate} because it exceeds {MAX_SIZE}")
                    continue
                variant_coordinate = variant_coordinate.as_external_explicit(genome_build)

                self._update_annotation(v, variant_coordinate, matcher, va_list, VariantAnnotation)
                self._update_annotation(v, variant_coordinate, matcher, vta_list, VariantTranscriptAnnotation)
                now = time.time()
                if now - last_update > 5:
                    last_update = now
                    print(f"{i} of {total} ({100 * i/total:.2f}%)")

            # Insert any remaining
            self._bulk_update(VariantAnnotation, va_list,1)
            self._bulk_update(VariantTranscriptAnnotation, vta_list,1)
