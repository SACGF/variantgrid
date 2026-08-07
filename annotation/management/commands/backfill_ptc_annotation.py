"""
Populates ptc_distance_codons / ptc_last_junction_distance / nmd_escape_status on rows
annotated before the PTC calculation existed (#579). hgvs_p and protein_position are already
stored, so the values are recoverable without reannotation.

Restricted to columns_version=5 VariantAnnotationVersions - that matches the
min_columns_version=5 gate on the column defs and where the inserter computes the values.

Only frameshift rows are touched. The inserter also writes NOT_APPLICABLE on everything else
(free, it's just another character in the COPY row), but writing it here would rewrite ~98% of
the partition for a value nothing reads - once backfilled_ptc is set, null carries the same
meaning to DamageNode.
"""
from django.core.management.base import BaseCommand

from annotation.models import VariantAnnotation, VariantAnnotationVersion
from annotation.models.models import VariantTranscriptAnnotation
from annotation.vcf_files.bulk_vep_vcf_annotation_inserter import (
    FRAMESHIFT_CONSEQUENCE,
    TranscriptGeometryCache,
    add_calculated_ptc,
)

BATCH_SIZE = 5000
PTC_FIELDS = ["ptc_distance_codons", "ptc_last_junction_distance", "nmd_escape_status"]


class Command(BaseCommand):
    def handle(self, *args, **options):
        for vav in VariantAnnotationVersion.objects.filter(columns_version=5):
            print(f"Updating {vav}...")
            # One cache per VAV - the geometry is per (transcript version, build)
            transcript_geometry_cache = TranscriptGeometryCache()
            for klass in (VariantTranscriptAnnotation, VariantAnnotation):
                qs = klass.objects.filter(version=vav, consequence__contains=FRAMESHIFT_CONSEQUENCE)
                num_frameshift = self._backfill_frameshifts(klass, qs, transcript_geometry_cache)
                print(f"  {klass.__name__}: {num_frameshift} frameshift rows")

            vav.backfilled_ptc = True
            vav.save(update_fields=["backfilled_ptc"])

    @staticmethod
    def _backfill_frameshifts(klass, qs, transcript_geometry_cache: TranscriptGeometryCache) -> int:
        qs = qs.select_related("variant__locus__ref", "variant__alt").order_by("pk")
        num_updated = 0
        batch = []
        for record in qs.iterator(chunk_size=BATCH_SIZE):
            ptc_data = {
                "consequence": record.consequence,
                "hgvs_p": record.hgvs_p,
                "protein_position": record.protein_position,
                "transcript_version_id": record.transcript_version_id,
            }
            ref = record.variant.locus.ref.seq
            alt = record.variant.alt.seq
            add_calculated_ptc(ptc_data, len(ref) - len(alt), transcript_geometry_cache)
            for field in PTC_FIELDS:
                setattr(record, field, ptc_data.get(field))
            batch.append(record)
            if len(batch) >= BATCH_SIZE:
                num_updated += _bulk_update(klass, batch)
                batch = []
        if batch:
            num_updated += _bulk_update(klass, batch)
        return num_updated


def _bulk_update(klass, batch: list) -> int:
    klass.objects.bulk_update(batch, PTC_FIELDS, batch_size=BATCH_SIZE)
    return len(batch)
