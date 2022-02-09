import numpy as np
from django.core.management.base import BaseCommand

from annotation.models import AnnotationRangeLock, VariantAnnotation, VariantTranscriptAnnotation
from snpdb.models import Variant, VariantZygosityCount, VariantCollectionRecord


class Command(BaseCommand):
    """
        Until 210622 - (PythonKnownVariantsImporter v.16) we used to insert a reference variant (alt='=') for each ALT
        We also didn't have a way to delete variants that are no longer referenced
    """
    def add_arguments(self, parser):
        # Usually an annotation range lock is 100k, so you'd expect 50k ref variants in there.
        # So steps=20 will look in a 5k range and probably delete 2.5k Variant and 2.5k VariantZygosityCount
        parser.add_argument('--steps', type=int, default=20, required=False,
                            help="Number of steps to take in between AnnotationRangeLock regions (which are ~100k)")
        parser.add_argument('--min-variant', type=int, required=False)

    def handle(self, *args, **options):
        # We want to do this in small batches - so use the variant annotation range locks which are all approx the same
        # size (even if a big gap between IDs)
        # Variants from different builds are mixed up together so just do 37 (won't be much of 38 that is not overlapping)
        # Any overlaps will be quick as not much to delete
        steps = options["steps"]
        min_variant = options["min_variant"]

        check_ref = True
        arl_qs = AnnotationRangeLock.objects.filter(version__genome_build__name='GRCh37')
        if min_variant:
            arl_qs = arl_qs.filter(min_variant__gt=min_variant)
        total = arl_qs.count()
        print(f"Deleting unused variants in {total} steps...")
        total_deleted = 0
        for i, range_lock in enumerate(arl_qs.order_by("max_variant")):
            perc = 100 * i / total
            print(f"{perc:.2f}% done - doing range lock {range_lock.pk}: {range_lock}")

            linspace = np.linspace(range_lock.min_variant_id, range_lock.max_variant_id, steps + 1).astype(int)
            for s in range(steps):
                start = linspace[s]
                end = linspace[s+1]
                num = end - start
                if num == 0:
                    continue
                print(f"{start} - {end} ({num})")
                # Skip the range lock min/max variant as that's protected so can't delete anyway
                # Also need to avoid deleting those that are within range in another build
                variants_in_range_qs = Variant.objects.filter(pk__gt=start, pk__lt=end)
                unused_variants_qs = variants_in_range_qs.filter(classification__isnull=True,
                                                                 clinvar__isnull=True,
                                                                 varianttag__isnull=True,
                                                                 variantallele__isnull=True,
                                                                 cohortgenotype__isnull=True,
                                                                 createdmanualvariant__isnull=True,
                                                                 min_variant__isnull=True,
                                                                 max_variant__isnull=True)
                unused_variant_ids = list(unused_variants_qs.values_list("pk", flat=True))
                if unused_variant_ids and check_ref:
                    print(f"{len(unused_variant_ids)=}")
                    num_ref = unused_variants_qs.filter(alt__seq='=').count()
                    print(f"check... {num_ref=}")
                    check_ref = False

                # If there is a VariantCollectionRecord, but no sample, the analysis won't work anyway
                for klass in [VariantCollectionRecord, VariantZygosityCount, VariantAnnotation, VariantTranscriptAnnotation]:
                    qs = klass.objects.filter(variant_id__gt=start, variant_id__lt=end, variant_id__in=unused_variant_ids)
                    num_deleted = qs._raw_delete(qs.db)
                    print(f"{klass}: deleted {num_deleted} records")

                variants_deleted = unused_variants_qs._raw_delete(unused_variants_qs.db)
                total_deleted += variants_deleted
                print(f"{variants_deleted=} ({total_deleted=})")

        print(f"Total deleted: {total_deleted}")
