import logging

import numpy as np
import pandas as pd
from django.conf import settings
from django.core.management.base import BaseCommand
from django.db import IntegrityError
from django.db.models import OuterRef, Subquery, F
from django.db.models.functions import Abs

from annotation.models import AnnotationRangeLock
from genes.hgvs import HGVSMatcher
from library.genomics.vcf_enums import VCFSymbolicAllele
from library.utils import md5sum_str
from snpdb.models import Variant, Locus, Sequence, GenomeBuild, VariantCoordinate


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
        parser.add_argument('--dry-run', action='store_true')
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
                                              length=len(symbolic_alt))  # This is wrong...
        return s

    @staticmethod
    def update_variants_in_range_make_symbolic(variant_qs, dry_run: bool):
        seq = {s: Sequence.objects.get(seq=s) for s in "GATC"}
        seq_del = Command._create_symbolic(VCFSymbolicAllele.DEL)
        seq_dup = Command._create_symbolic(VCFSymbolicAllele.DUP)
        changed_rows = []

        # Make sure we don't re-do stuff
        variant_qs = variant_qs.exclude(alt__in=[seq_del, seq_dup])

        # Change to dels
        for v in variant_qs.filter(locus__ref__length__gte=1000, alt__length=1).select_related("locus__ref", "alt"):
            ref = v.locus.ref.seq[0]
            alt = v.alt.seq[0]
            if ref != alt:
                continue  # delins

            changed_data = {
                "variant_id": v.pk,
                "contig": v.locus.contig.name,
                "position": v.locus.position,
                "old_ref": v.locus.ref.seq,
                "old_alt": v.alt.seq,
                "genome_build": next(iter(v.genome_builds)).name,
            }

            if not dry_run:
                locus, created = Locus.objects.get_or_create(contig=v.locus.contig, position=v.locus.position, ref=seq[ref])
                try:
                    v.locus = locus
                    v.alt = seq_del
                    v.save()
                except IntegrityError as e:
                    print(e)
                    existing = Variant.objects.get(locus=locus, end=v.end, alt=seq_del)
                    msg = f"Variant pk={v.pk} is a dupe of pk={existing.pk}"
                    logging.error(msg)
                    continue

            changed_data["end"] = v.end
            changed_data["new_ref"] = ref
            changed_data["new_alt"] = seq_del.seq
            changed_rows.append(changed_data)

        new_dups = []

        # DUPS - need to check if an actual dup (not just an ins) - easiest to do via HGVS
        # Need to retrieve genome sequence so have to handle builds separately
        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build)

            build_qs = variant_qs.filter(Variant.get_contigs_q(genome_build))
            build_qs = build_qs.filter(locus__ref__length=1, alt__length__gte=settings.VARIANT_SYMBOLIC_ALT_SIZE)
            for v in build_qs.select_related("locus__ref", "alt"):
                ref = v.locus.ref.seq[0]
                alt = v.alt.seq[0]
                if ref != alt:
                    continue  # quick ins reject

                try:
                    hgvs_variant = matcher.variant_coordinate_to_hgvs_variant(v.coordinate)
                    if hgvs_variant.mutation_type == 'dup':
                        changed_data = {
                            "variant_id": v.pk,
                            "contig": v.locus.contig.name,
                            "position": v.locus.position,
                            "old_ref": v.locus.ref.seq,
                            "old_alt": v.alt.seq,
                            "genome_build": next(iter(v.genome_builds)).name,
                            "end": v.end,
                            "new_ref": v.locus.ref.seq,
                            "new_alt": seq_dup.seq,
                        }
                        changed_rows.append(changed_data)

                        # Will do this in bulk
                        v.alt_id = seq_dup.pk
                        new_dups.append(v)
                except Exception as e:
                    logging.error("Couldn't convert %s to HGVS: %s", v, e)

        if not dry_run and len(new_dups):
            Variant.objects.bulk_update(new_dups, ["alt_id"], batch_size=2000)

        return changed_rows

    def handle(self, *args, **options):
        # We want to do this in small batches - so use the variant annotation range locks which are all approx the same
        # size (even if a big gap between IDs)
        # Variants from different builds are mixed up together - we just want the biggest one
        steps = options["steps"]
        replace = options["replace"]
        dry_run = options["dry_run"]
        min_variant = options["min_variant"]

        changed_rows = []
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
                changed_rows.extend(self.update_variants_in_range_make_symbolic(variant_qs, dry_run))
                last_max = end

            if changed_rows:
                build_matchers = {
                    "GRCh37": HGVSMatcher(GenomeBuild.grch37()),
                    "GRCh38": HGVSMatcher(GenomeBuild.grch38()),
                }

                different = 0
                # calculate
                for row in changed_rows:
                    vc_old = VariantCoordinate.from_explicit_no_svlen(chrom=row["contig"], start=row["position"],
                                                                      ref=row["old_ref"], alt=row["old_alt"])
                    vc_new = VariantCoordinate(chrom=row["contig"], start=row["position"],
                                               end=row["end"], ref=row["new_ref"], alt=row["new_alt"])

                    matcher: HGVSMatcher = build_matchers[row["genome_build"]]
                    try:
                        g_hgvs_old = matcher.variant_coordinate_to_g_hgvs(vc_old)
                    except Exception as e:
                        g_hgvs_old = "ERROR: " + str(e)

                    try:
                        g_hgvs_new = matcher.variant_coordinate_to_g_hgvs(vc_new)
                    except Exception as e:
                        g_hgvs_new = "ERROR: " + str(e)

                    row["g_hgvs_old"] = g_hgvs_old
                    row["g_hgvs_new"] = g_hgvs_new
                    same = g_hgvs_new == g_hgvs_old
                    row["same"] = same
                    if not same:
                        different += 1

                if different:
                    print(f"There were {different} differences!!")
                df = pd.DataFrame.from_records(changed_rows)
                df.to_csv("variant_symbolic_alt_changes.csv", index=False)

        # Delete any orphaned locus/sequences
        if not dry_run:
            Locus.objects.filter(variant__isnull=True).delete()
            Sequence.objects.filter(locus__isnull=True, variant__isnull=True).delete()
