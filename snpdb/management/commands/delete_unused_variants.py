from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import Max

from annotation.models import VariantAnnotation, VariantTranscriptAnnotation, VariantGeneOverlap
from snpdb.models import Variant, VariantZygosityCount, VariantCollectionRecord
from upload.models import ModifiedImportedVariant


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

            # Any model with a FK to Variant must be listed here (to keep the variant) or in the
            # raw-delete loop below (to delete the dependent record) - otherwise the variant delete
            # below will fail with an IntegrityError as _raw_delete bypasses Django's on_delete.
            unused_variants_qs = Variant.objects.filter(pk__gte=start, pk__lte=end,
                                                        classification__isnull=True,
                                                        clinvar__isnull=True,
                                                        varianttag__isnull=True,
                                                        variantallele__isnull=True,
                                                        cohortgenotype__isnull=True,
                                                        commonvariantclassified__isnull=True,
                                                        createdmanualvariant__isnull=True,
                                                        candidate__isnull=True,
                                                        nodevariant__isnull=True,
                                                        intersectionnode__isnull=True,
                                                        importedalleleinfo__isnull=True,
                                                        resolvedvariantinfo__isnull=True,
                                                        variantwiki__isnull=True,
                                                        uploadedvcf__isnull=True,
                                                        min_variant__isnull=True,
                                                        max_variant__isnull=True,
                                                        allvariantsnode__isnull=True)
            unused_variant_ids = list(unused_variants_qs.values_list("pk", flat=True))
            if not unused_variant_ids:
                continue
            if check_ref:
                print(f"{len(unused_variant_ids)=}")
                num_ref = unused_variants_qs.filter(alt__seq='=').count()
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
