import os

import numpy as np
from django.conf import settings
from django.core.management.base import BaseCommand
from django.db.models import OuterRef, Subquery, F
from django.db.models.functions import Length, Abs

from annotation.models import AnnotationRangeLock
from library.utils import mk_path
from snpdb.models import Variant


class Command(BaseCommand):
    """
        This is to fix variant end after we changed to using SVLEN

        @see https://github.com/SACGF/variantgrid/issues/990
    """
    def add_arguments(self, parser):
        # Usually an annotation range lock is 100k, so you'd expect 50k ref variants in there.
        # So steps=20 will look in a 5k range
        parser.add_argument('--steps', type=int, default=20, required=False,
                            help="Number of steps to take in between AnnotationRangeLock regions (which are ~100k)")
        parser.add_argument('--replace', action='store_true')
        parser.add_argument('--min-variant', type=int, required=False)

    @staticmethod
    def update_variants_in_range_fix_end(variant_qs):
        # These are all non-symbolic, so can just use ref length
        non_symbolic_variant_qs = variant_qs.filter(svlen__isnull=True)
        calc_end = F("locus__position") + Length("locus__ref__seq") - 1
        non_symbolic_variant_subquery = Variant.objects.filter(pk=OuterRef("pk")).annotate(calc_end=calc_end).values("calc_end")[:1]
        non_symbolic_variant_qs.update(end=Subquery(non_symbolic_variant_subquery))

        symbolic_variant_qs = variant_qs.filter(svlen__isnull=False)
        calc_end = F("locus__position") + Abs("svlen")
        symbolic_variant_subquery = Variant.objects.filter(pk=OuterRef("pk")).annotate(calc_end=calc_end).values("calc_end")[:1]
        symbolic_variant_qs.update(end=Subquery(symbolic_variant_subquery))

    def handle(self, *args, **options):
        # We want to do this in small batches - so use the variant annotation range locks which are all approx the same
        # size (even if a big gap between IDs)
        # Variants from different builds are mixed up together - we just want the biggest one
        steps = options["steps"]
        min_variant = options["min_variant"]

        # This can take a few days, so we'll write the variant as we go, so we can resume without any troubles
        migrations_dir = os.path.join(settings.PRIVATE_DATA_ROOT, "migrations")
        mk_path(migrations_dir)
        progress_file = os.path.join(migrations_dir, "one_off_fix_variant_end_progress_v2.txt")

        highest_av = AnnotationRangeLock.objects.order_by("-max_variant").first()
        arl_qs = AnnotationRangeLock.objects.filter(version=highest_av.version)

        if min_variant is None:
            # Try loading from CSV file
            if os.path.exists(progress_file):
                with open(progress_file) as f:
                    line = f.read()
                    min_variant = int(line)
                    print(f"Starting from {min_variant=} read from progress file...")

        last_max = 0  # Use previous max to range lock max as it may have skipped some as it only has 1 build
        if min_variant:
            arl_qs = arl_qs.filter(min_variant__gt=min_variant)
            last_max = min_variant

        total = arl_qs.count()
        print(f"Adding start/end to variants in {total} steps...")
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

                variant_qs = Variant.objects.filter(pk__gte=start, pk__lte=end)
                self.update_variants_in_range_fix_end(variant_qs)

                last_max = end

            # Record progress - we'll redo the last one just to be sure
            with open(progress_file, "w") as f:
                f.write(str(range_lock.min_variant_id))
