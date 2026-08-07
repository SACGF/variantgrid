from collections import Counter
from functools import cache

from django.core.management.base import BaseCommand
from django.db.models import Q
from django.db.models.functions import Length

from annotation.models import AnnotationRangeLock, ClinVar
from genes.hgvs import HGVSMatcher
from snpdb.models import GenomeBuild, Locus, Sequence, Variant


@cache
def _get_sequence(seq: str) -> Sequence:
    """ Sequences are created on first sight, so eg <DEL> won't exist yet on a legacy DB """
    sequence, _ = Sequence.objects.get_or_create(seq=seq)
    return sequence


class Command(BaseCommand):
    """
        There was a bug where SVs were not converted to symbolic for a while in #982

        This will try and fix them
    """
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        long_sequences = Sequence.objects.all().annotate(seq_length=Length("seq")).filter(seq_length__gte=1000)
        # select_related as v.coordinate reads locus/contig/ref/alt
        long_variants = Variant.objects.filter(Q(locus__ref__in=long_sequences) | Q(alt__in=long_sequences)) \
            .select_related("locus__contig", "locus__ref", "alt")
        print(f"Long variant count = {long_variants.count()}")

        dry_run = options["dry_run"]

        for genome_build in GenomeBuild.builds_with_annotation():
            print(f"Genome build {genome_build}")
            self._find_bad_symbolic_via_clinvar(dry_run, genome_build, "<DEL>")  # DELISN stored as DEL
            self._find_bad_symbolic_via_clinvar(dry_run, genome_build, "<DUP>")  # INS stored as DUP

            stats = Counter()
            matcher = HGVSMatcher.instance(genome_build)
            q_contig = Variant.get_contigs_q(genome_build)
            for v in long_variants.filter(q_contig):
                vc = v.coordinate.as_internal_symbolic(genome_build)

                try:
                    old_hgvs = matcher.variant_coordinate_to_g_hgvs(v.coordinate)
                    new_hgvs = matcher.variant_coordinate_to_g_hgvs(vc)
                except Exception as e:
                    stats[f"g.HGVS check skipped ({type(e).__name__})"] += 1
                else:
                    if old_hgvs != new_hgvs:
                        stats["g.HGVS mismatch - old != new"] += 1

                if not vc.is_symbolic:
                    stats["left explicit - no symbolic form"] += 1
                    continue

                try:
                    existing = Variant.get_from_variant_coordinate(vc, genome_build)
                    if v == existing:
                        ref_length = len(v.locus.ref)
                        if ref_length > 1:
                            # It's already there as existing - but didn't have the ref trimmed down
                            stats["ref trimmed on existing symbolic"] += 1
                            if not dry_run:
                                new_ref = _get_sequence(v.locus.ref.seq[0])
                                v.locus = Locus.objects.get_or_create(contig=v.locus.contig, position=vc.position,
                                                                      ref=new_ref)[0]
                                v.end = vc.end
                                v.save()
                        else:
                            stats["already symbolic - no change"] += 1
                    else:
                        # I think these came in via ClinVar import - so will try to get rid of it, if it isn't referenced
                        # by anything else
                        stats["merged into existing symbolic"] += 1
                        if not dry_run:
                            self._merge_variant_dupe(v, existing, stats)
                except Variant.DoesNotExist:
                    if len(vc.ref) > 1:
                        raise ValueError(f"{v.pk} had ref length of {len(vc.ref)}")

                    stats[f"converted to {vc.alt}"] += 1
                    if not dry_run:
                        new_ref = _get_sequence(vc.ref)
                        v.locus = Locus.objects.get_or_create(contig=v.locus.contig, position=vc.position, ref=new_ref)[0]
                        v.alt = _get_sequence(vc.alt)
                        v.svlen = vc.svlen
                        v.end = vc.end  # Stored calculated field - nothing recalcs it on save
                        v.save()

            self._print_stats(f"{genome_build} long variants", stats)

        # This will take ages on some systems...
        # Locus.objects.filter(variant__isnull=True).delete()

    @staticmethod
    def _print_stats(label: str, stats: Counter):
        print(f"{label}:")
        if not stats:
            print("  (nothing)")
        for reason, count in stats.most_common():
            print(f"  {count}\t{reason}")

    def _merge_variant_dupe(self, dupe_variant, original_variant, stats: Counter) -> int:
        if dupe_variant.cohortgenotype_set.exists():
            stats["NOT merged - has cohort genotype data"] += 1
            return 0

        dupe_variant.clinvar_set.all().update(variant=original_variant)
        AnnotationRangeLock.release_variant(dupe_variant)
        # Classification/Variant tags protects Variant FK, so won't delete if that exists
        try:
            dupe_variant.delete()
            return 1
        except Exception as e:
            stats[f"NOT merged - protected ({type(e).__name__})"] += 1
        return 0

    def _find_bad_symbolic_via_clinvar(self, dry_run: bool, genome_build, alt_seq):
        """ Due to VariantCoordinate.as_symbolic_variant bug - some historical data was incorrectly imported """

        stats = Counter()

        clinvar_variation_del = list(
            ClinVar.objects.filter(version__genome_build=genome_build, variant__alt__seq=alt_seq).values_list(
                "clinvar_variation_id", flat=True))

        clinvar_variation_original = {}

        for cv in ClinVar.objects.filter(version__genome_build=genome_build,
                                         clinvar_variation_id__in=clinvar_variation_del).exclude(
                variant__alt__seq=alt_seq):
            clinvar_variation_original[cv.clinvar_variation_id] = cv.variant

        clinvar_variation_bad = {}
        for cv in ClinVar.objects.filter(version__genome_build=genome_build,
                                         clinvar_variation_id__in=clinvar_variation_original,
                                         variant__alt__seq=alt_seq):
            clinvar_variation_bad[cv.clinvar_variation_id] = cv.variant

        if not clinvar_variation_bad:
            return

        for clinvar_variation_id, bad_variant in clinvar_variation_bad.items():
            original_variant = clinvar_variation_original[clinvar_variation_id]
            stats["merged into original representation"] += 1
            if not dry_run:
                self._merge_variant_dupe(bad_variant, original_variant, stats)

        self._print_stats(f"{genome_build} clinvar {alt_seq} ({len(clinvar_variation_del)} records, "
                          f"{len(clinvar_variation_bad)} bad)", stats)
