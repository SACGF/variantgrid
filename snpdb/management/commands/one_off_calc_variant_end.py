import logging
from collections import Counter

import numpy as np
from django.core.management.base import BaseCommand
from django.db.models import OuterRef, Subquery, F
from django.db.models.functions import Abs

from annotation.models import AnnotationRangeLock
from genes.hgvs import HGVSMatcher
from library.utils import md5sum_str
from snpdb.models import Variant, Locus, Sequence, GenomeBuild


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
        parser.add_argument('--replace', action='store_true')
        parser.add_argument('--min-variant', type=int, required=False)

    @staticmethod
    def update_variants_in_range_add_end(variant_qs, replace=False):
        # Set End
        calc_end = F("locus__position") + Abs(F("locus__ref__length") - F("alt__length"))
        variant_subquery = Variant.objects.filter(pk=OuterRef("pk")).annotate(calc_end=calc_end).values("calc_end")[:1]
        if not replace:
            variant_qs = variant_qs.filter(end__isnull=True)
        variant_qs.update(end=Subquery(variant_subquery))

    @staticmethod
    def _create_symbolic(symbolic_alt):
        s, _ = Sequence.objects.get_or_create(seq=symbolic_alt,
                                              seq_md5_hash=md5sum_str(symbolic_alt),
                                              length=len(symbolic_alt)) # This is wrong...
        return s

    @staticmethod
    def update_variants_in_range_make_symbolic(variant_qs):
        seq = {s: Sequence.objects.get(seq=s) for s in "GATC"}
        seq_del = Command._create_symbolic("<DEL>")
        seq_dup = Command._create_symbolic("<DUP>")
        changed = Counter()

        # Make sure we don't re-do stuff
        variant_qs = variant_qs.exclude(alt__in=[seq_del, seq_dup])

        # Change to dels
        for v in variant_qs.filter(locus__ref__length__gte=1000, alt__length=1).select_related("locus__ref", "alt"):
            ref = v.locus.ref.seq[0]
            alt = v.alt.seq[0]
            if ref != alt:
                continue  # delins

            locus, created = Locus.objects.get_or_create(contig=v.locus.contig, position=v.locus.position, ref=seq[ref])
            v.locus = locus
            v.alt = seq_del
            v.save()
            changed["<DEL>"] += 1

        new_dups = []

        # DUPS - need to check if an actual dup (not just an ins) - easiest to do via HGVS
        # Need to retrieve genome sequence so have to handle builds separately
        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build)

            build_qs = variant_qs.filter(Variant.get_contigs_q(genome_build))
            for v in build_qs.filter(locus__ref__length=1, alt__length__gte=1000).select_related("locus__ref", "alt"):
                ref = v.locus.ref.seq[0]
                alt = v.alt.seq[0]
                if ref != alt:
                    continue  # quick ins reject

                try:
                    hgvs_variant = matcher.variant_coordinate_to_hgvs_variant(v.coordinate)
                    if hgvs_variant.mutation_type == 'dup':
                        # Will do this in bulk
                        v.alt_id = seq_dup.pk
                        new_dups.append(v)
                except Exception as e:
                    logging.error("Couldn't convert %s to HGVS: %s", v, e)

        if len(new_dups):
            changed["<DUP>"] += len(new_dups)
            Variant.objects.bulk_update(new_dups, ["alt_id"], batch_size=2000)

        return changed

    def handle(self, *args, **options):
        # We want to do this in small batches - so use the variant annotation range locks which are all approx the same
        # size (even if a big gap between IDs)
        # Variants from different builds are mixed up together - we just want the biggest one
        steps = options["steps"]
        replace = options["replace"]
        min_variant = options["min_variant"]

        changed = Counter()
        missing_end = Variant.objects.filter(end__isnull=True).exists()
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
                variant_qs = Variant.objects.filter(pk__gte=start, pk__lte=end)
                if missing_end or replace:
                    self.update_variants_in_range_add_end(variant_qs, replace=replace)
                changed.update(self.update_variants_in_range_make_symbolic(variant_qs))
                last_max = end

            if changed:
                print(f"Changed: {changed}")

        # Now, go through and delete any orphaned:
        # * Locus w/o Variant
        # * Sequence w/o Variant or Locus

