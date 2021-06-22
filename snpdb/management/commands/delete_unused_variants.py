from django.core.management.base import BaseCommand
import logging

from annotation.models import VariantAnnotationVersion, AnnotationRangeLock
from snpdb.models import GenomeBuild, Variant


class Command(BaseCommand):
    """
        Until 210622 - (PythonKnownVariantsImporter v.16) we used to insert a reference variant (alt='=') for each ALT
        We also didn't have a way to delete variants that are no longer referenced
    """
    def handle(self, *args, **options):
        # We want to do this in small batches - so use the variant annotation range locks which are all approx the same
        # size (even if a big gap between IDs)
        # Variants from different builds are mixed up together so just do all builds
        # Any overlaps will be quick as not much to delete

        arl_qs = AnnotationRangeLock.objects.all()
        total = arl_qs.count()
        print(f"Deleting unused variants in {total} steps...")
        total_deleted = 0
        for i, range_lock in enumerate(arl_qs.order_by("max_variant")):
            perc = 100 * i / total
            print(f"{perc:.2f}% done - doing step: {range_lock}")

            # Skip the range lock min/max variant as that's protected so can't delete anyway
            qs = Variant.objects.filter(pk__gt=range_lock.min_variant_id, pk__lt=range_lock.max_variant_id)
            qs = qs.filter(classification__isnull=True,
                           clinvar__isnull=True,
                           varianttag__isnull=True,
                           variantallele__isnull=True,
                           cohortgenotype__isnull=True,
                           # Avoid deleting those that are are within range in another build
                           min_variant__isnull=True,
                           max_variant__isnull=True)
            _, details = qs.delete()
            num_deleted = details.get('snpdb.Variant', 0)
            print(details)
            total_deleted += num_deleted

        print(f"Total deleted: {total_deleted}")
        print("You should probably run:")
        print("python3.8 manage.py load_variants_hash_in_redis --clear")
        # Convert any ref=alt sequences into alt='='

