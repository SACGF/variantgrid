from django.core.management.base import BaseCommand
import logging
import numpy as np

from annotation.models import VariantAnnotationVersion, AnnotationRangeLock, VariantAnnotation
from snpdb.models import GenomeBuild, Variant, VariantZygosityCount


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

    def handle(self, *args, **options):
        # We want to do this in small batches - so use the variant annotation range locks which are all approx the same
        # size (even if a big gap between IDs)
        # Variants from different builds are mixed up together so just do all builds
        # Any overlaps will be quick as not much to delete
        steps = options["steps"]

        arl_qs = AnnotationRangeLock.objects.all()
        total = arl_qs.count()
        print(f"Deleting unused variants in {total} steps...")
        total_deleted = 0
        for i, range_lock in enumerate(arl_qs.order_by("max_variant")):
            perc = 100 * i / total
            print(f"{perc:.2f}% done - doing step: {range_lock}")

            linspace = np.linspace(range_lock.min_variant_id, range_lock.max_variant_id, steps + 1).astype(int)
            for s in range(steps):
                start = linspace[s]
                end = linspace[s+1]
                print(f"{start} - {end} ({end-start})")
                # Skip the range lock min/max variant as that's protected so can't delete anyway
                # Also need to avoid deleting those that are are within range in another build
                variants_in_range_qs = Variant.objects.filter(pk__gt=start, pk__lt=end)
                unused_variants_qs = variants_in_range_qs.filter(classification__isnull=True,
                                                                 clinvar__isnull=True,
                                                                 varianttag__isnull=True,
                                                                 variantallele__isnull=True,
                                                                 cohortgenotype__isnull=True,
                                                                 variantcollectionrecord__isnull=True,
                                                                 min_variant__isnull=True,
                                                                 max_variant__isnull=True)
                variants_deleted = unused_variants_qs._raw_delete(unused_variants_qs.db)
                print(f"{variants_deleted=}")
                vzc_qs = VariantZygosityCount.objects.filter(variant_id__gt=start, variant_id__lt=end)
                vzc_qs = vzc_qs.exclude(variant__in=variants_in_range_qs)
                zygosity_count_deleted = vzc_qs._raw_delete(vzc_qs.db)
                print(f"{zygosity_count_deleted=}")

                va_qs = VariantAnnotation.objects.filter(variant_id__gt=start, variant_id__lt=end)
                va_qs = va_qs.exclude(variant__in=variants_in_range_qs)
                annotation_deleted = va_qs._raw_delete(va_qs.db)
                print(f"{annotation_deleted=}")
                total_deleted += variants_deleted

        print(f"Total deleted: {total_deleted}")
        print("You should probably run:")
        print("python3.8 manage.py load_variants_hash_in_redis --clear")
        # Convert any ref=alt sequences into alt='='

