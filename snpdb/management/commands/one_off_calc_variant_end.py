import numpy as np
from django.core.management.base import BaseCommand
from django.db.models import OuterRef, Subquery
from django.db.models.functions import Greatest

from annotation.models import AnnotationRangeLock
from snpdb.models import Variant


class Command(BaseCommand):
    """
        We added the Variant end int fields for CNV work, leaving them nullable
        This command populates it on legacy imports so all Variant objects are eventually populated
    """
    def add_arguments(self, parser):
        # Usually an annotation range lock is 100k, so you'd expect 50k ref variants in there.
        # So steps=20 will look in a 5k range
        parser.add_argument('--steps', type=int, default=20, required=False,
                            help="Number of steps to take in between AnnotationRangeLock regions (which are ~100k)")
        parser.add_argument('--min-variant', type=int, required=False)

    @staticmethod
    def update_variants_in_range(start, end):
        variant_subquery = Variant.objects.filter(pk=OuterRef("pk")) \
                               .annotate(longest=Greatest("locus__ref__length", "alt__length")) \
                               .values_list("longest")[:1]
        Variant.objects.filter(pk__gte=start, pk__lte=end) \
            .filter(end__isnull=True).update(end=Subquery(variant_subquery))

    def handle(self, *args, **options):
        # We want to do this in small batches - so use the variant annotation range locks which are all approx the same
        # size (even if a big gap between IDs)
        # Variants from different builds are mixed up together - we just want the biggest one
        steps = options["steps"]
        min_variant = options["min_variant"]

        highest_av = AnnotationRangeLock.objects.order_by("-max_variant").first()
        arl_qs = AnnotationRangeLock.objects.filter(version=highest_av.version)
        if min_variant:
            arl_qs = arl_qs.filter(min_variant__gt=min_variant)
        total = arl_qs.count()
        print(f"Adding start/end to variants in {total} steps...")
        last_max = 0  # Use previous max to range lock max as it may have skipped some as it only has 1 build
        for i, range_lock in enumerate(arl_qs.order_by("max_variant")):
            perc = 100 * i / total
            print(f"{perc:.2f}% done - doing range lock {range_lock.pk}: {last_max} - {range_lock.max_variant_id}")

            linspace = np.linspace(last_max, range_lock.max_variant_id, steps + 1).astype(int)
            for s in range(steps):
                start = linspace[s]
                end = linspace[s+1]
                num = end - start
                if num == 0:
                    continue
                self.update_variants_in_range(start, end)
                last_max = end
