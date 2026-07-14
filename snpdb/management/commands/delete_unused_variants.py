from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import Max

from analysis.models import AllVariantsNode, Candidate, IntersectionNode, NodeVariant, VariantTag
from annotation.models import VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap, \
    ClinVar, CreatedManualVariant, AnnotationRangeLock
from classification.models import Classification, ImportedAlleleInfo, ResolvedVariantInfo
from snpdb.models import Variant, VariantZygosityCount, VariantCollectionRecord, \
    CohortGenotype, CommonVariantClassified, VariantAllele, VariantWiki
from upload.models import ModifiedImportedVariant, UploadedVCFPipelineMaxVariant


# Every model with a FK to Variant whose presence means the variant is still referenced and must be
# kept, as (model, fk_column). A variant is "unused" only if no row in any of these (or in
# PRELOAD_VARIANT_RELATIONS below) points at it. Keep this in sync with FKs to Variant - a missing
# entry would delete still-referenced variants; the raw-delete loop in handle() covers the
# derived/regenerable models we instead delete along with the variant.
# These are the large (≈one row per variant) tables, range-scanned once per batch.
KEEP_VARIANT_RELATIONS = [
    (Classification, "variant_id"),
    (ClinVar, "variant_id"),
    (VariantTag, "variant_id"),
    (VariantAllele, "variant_id"),
    (CohortGenotype, "variant_id"),
    (CommonVariantClassified, "variant_id"),
    (CreatedManualVariant, "variant_id"),
    (Candidate, "variant_id"),
    (NodeVariant, "variant_id"),
    (ImportedAlleleInfo, "matched_variant_id"),
    (ResolvedVariantInfo, "variant_id"),
    (VariantWiki, "variant_id"),
]

# Small tables - a bounded number of rows (≈one per annotation range / node / uploaded file, not per
# variant) - whose entire set of referenced variant PKs we load once up front, rather than
# range-scanning every batch.
PRELOAD_VARIANT_RELATIONS = [
    (AnnotationRangeLock, "min_variant_id"),
    (AnnotationRangeLock, "max_variant_id"),
    (AllVariantsNode, "max_variant_id"),
    (IntersectionNode, "hgvs_variant_id"),
    (UploadedVCFPipelineMaxVariant, "max_variant_id"),
]


class Command(BaseCommand):
    """
        Until 210622 - (PythonKnownVariantsImporter v.16) we used to insert a reference variant (alt='=') for each ALT
        We also didn't have a way to delete variants that are no longer referenced
    """
    def add_arguments(self, parser):
        parser.add_argument('--batch-size', type=int, default=5000, required=False,
                            help="Number of (actual) variants to examine per step")
        parser.add_argument('--min-variant', type=int, required=False,
                            help="Resume point - only examine variants with pk greater than this")

    def handle(self, *args, **options):
        # Walk the real Variant PKs with a keyset cursor (pk > last_pk ORDER BY pk LIMIT batch_size). This is an
        # index-only scan over the PK B-tree, so the (often huge) gaps between PKs are skipped for free - every
        # step does exactly batch_size variants of work. Variants from all genome builds share the one PK space,
        # so a single pass covers every build (no per-build range-lock iteration needed).
        batch_size = options["batch_size"]
        last_pk = options["min_variant"] or 0
        max_pk = Variant.objects.aggregate(m=Max("pk"))["m"] or 0

        check_ref = True
        total_deleted = 0

        # Load the small tables' referenced variant PKs once up front (see PRELOAD_VARIANT_RELATIONS)
        # rather than range-scanning them every batch.
        preloaded_referenced_ids = set()
        for klass, fk in PRELOAD_VARIANT_RELATIONS:
            preloaded_referenced_ids.update(klass.objects.values_list(fk, flat=True))

        while True:
            batch_pks = list(Variant.objects.filter(pk__gt=last_pk).order_by("pk")
                             .values_list("pk", flat=True)[:batch_size])
            if not batch_pks:
                break
            # batch_pks are the first batch_size ordered PKs, so nothing exists strictly between
            # start and end that isn't in the batch - the range selects exactly this batch.
            start, end = batch_pks[0], batch_pks[-1]
            last_pk = end
            perc = 100 * end / max_pk if max_pk else 100
            print(f"{perc:.2f}% done - variants {start} - {end} ({len(batch_pks)})")

            # A variant is unused if nothing references it - neither a preloaded small-table PK nor a
            # row in any KEEP_VARIANT_RELATIONS table. Rather than one query LEFT JOINing all ~18
            # tables (which blows past Postgres's join_collapse_limit of 8 and degrades into
            # full-table hash joins every batch), probe each big table with a single indexed range
            # scan over [start, end] and subtract the referenced PKs in Python. Each scan returns at
            # most ~batch_size ids, so this stays cheap.
            # Any model with a FK to Variant must be in KEEP_VARIANT_RELATIONS / PRELOAD_VARIANT_RELATIONS
            # (to keep the variant) or in the raw-delete loop below (to delete the dependent record) -
            # otherwise the variant delete will fail with an IntegrityError as _raw_delete bypasses
            # Django's on_delete.
            referenced_ids = set(preloaded_referenced_ids)
            for klass, fk in KEEP_VARIANT_RELATIONS:
                qs = klass.objects.filter(**{f"{fk}__gte": start, f"{fk}__lte": end})
                referenced_ids.update(qs.values_list(fk, flat=True))
            unused_variant_ids = [pk for pk in batch_pks if pk not in referenced_ids]
            if not unused_variant_ids:
                continue
            if check_ref:
                print(f"{len(unused_variant_ids)=}")
                num_ref = Variant.objects.filter(pk__in=unused_variant_ids, alt__seq='=').count()
                print(f"check... {num_ref=}")
                check_ref = False

            with transaction.atomic():
                # If there is a VariantCollectionRecord, but no sample, the analysis won't work anyway
                # VariantGeneOverlap/annotation rows are derived data that can be regenerated.
                for klass in [VariantCollectionRecord, VariantZygosityCount,
                              VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap,
                              ModifiedImportedVariant]:
                    qs = klass.objects.filter(variant_id__in=unused_variant_ids)
                    num_deleted = qs._raw_delete(qs.db)
                    print(f"{klass}: deleted {num_deleted} records")

                # Reuse the already-computed PKs rather than re-running the expensive unused filter
                delete_qs = Variant.objects.filter(pk__in=unused_variant_ids)
                variants_deleted = delete_qs._raw_delete(delete_qs.db)
            total_deleted += variants_deleted
            print(f"{variants_deleted=} ({total_deleted=})")

        print(f"Total deleted: {total_deleted}")
