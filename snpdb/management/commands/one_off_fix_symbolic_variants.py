from django.core.management.base import BaseCommand
from django.db.models import Q
from django.db.models.functions import Length

from annotation.models import AnnotationRangeLock, ClinVar
from genes.hgvs import HGVSMatcher
from snpdb.models import Variant, Sequence, GenomeBuild, Locus


class Command(BaseCommand):
    """
        There was a bug where SVs were not converted to symbolic for a while in #982

        This will try and fix them
    """
    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true')

    def handle(self, *args, **options):
        long_sequences = Sequence.objects.all().annotate(seq_length=Length("seq")).filter(seq_length__gte=1000)
        long_variants = Variant.objects.filter(Q(locus__ref__in=long_sequences) | Q(alt__in=long_sequences))

        base_lookup = {s: Sequence.objects.get(seq=s) for s in ["G", "A", "T", "C", "<DEL>", "<DUP>"]}
        not_symbolic = []
        dry_run = options["dry_run"]

        for genome_build in GenomeBuild.builds_with_annotation():
            self._find_bad_symbolic_via_clinvar(dry_run, genome_build, "<DEL>")  # DELISN stored as DEL
            self._find_bad_symbolic_via_clinvar(dry_run, genome_build, "<DUP>")  # INS stored as DUP
            num_deleted = 0

            matcher = HGVSMatcher(genome_build)
            q_contig = Variant.get_contigs_q(genome_build)
            for v in long_variants.filter(q_contig):
                vc = v.coordinate.as_internal_symbolic(genome_build)

                try:
                    old_hgvs = matcher.variant_coordinate_to_g_hgvs(v.coordinate)
                    new_hgvs = matcher.variant_coordinate_to_g_hgvs(vc)

                    if old_hgvs != new_hgvs:
                        raise ValueError(f"Variant: {v.pk} old: {old_hgvs} != {new_hgvs}")
                except Exception as e:
                    print(f"Could not test {v.pk} g.hgvs: {e}")

                if not vc.is_symbolic():
                    not_symbolic.append(v)
                    continue

                try:
                    existing = Variant.get_from_variant_coordinate(vc, genome_build)
                    if v == existing:
                        ref_length = len(v.locus.ref)
                        if ref_length > 1:
                            # It's already there as existing - but didn't have the ref trimmed down
                            print(f"Fixing {v.pk} ({v.alt}) had ref sequence length of {ref_length}")
                            if not dry_run:
                                new_ref = base_lookup[v.locus.ref.seq[0]]
                                v.locus = Locus.objects.get_or_create(contig=v.locus.contig, position=vc.position,
                                                                      ref=new_ref)[0]
                                v.save()
                    else:
                        print(f"Deleting Variant={v.pk}")
                        # I think these came in via ClinVar import - so will try to get rid of it, if it isn't referenced
                        # by anything else
                        if not dry_run:
                            num_deleted += self._merge_variant_dupe(v, existing)
                except Variant.DoesNotExist:
                    print(f"Fixing {v.pk}: {v} - changing {vc.alt}")
                    if len(vc.ref) > 1:
                        raise ValueError(f"{v.pk} had ref length of {len(vc.ref)}")

                    if not dry_run:
                        new_ref = base_lookup[vc.ref]
                        v.locus = Locus.objects.get_or_create(contig=v.locus.contig, position=vc.position, ref=new_ref)[0]
                        v.alt = base_lookup[vc.alt]
                        v.svlen = vc.svlen
                        v.save()

            print(f"{GenomeBuild} Merged/Deleted {num_deleted} variants")

        print(f"Non-symbolic not converted: {len(not_symbolic)}")
        # This will take ages on some systems...
        # Locus.objects.filter(variant__isnull=True).delete()

    def _merge_variant_dupe(self, dupe_variant, original_variant) -> int:
        if dupe_variant.cohortgenotype_set.exists():
            print(f"Not deleting {dupe_variant.pk} as it has cohort genotype data")
            return

        dupe_variant.clinvar_set.all().update(variant=original_variant)
        AnnotationRangeLock.release_variant(dupe_variant)
        # Classification/Variant tags protects Variant FK, so won't delete if that exists
        try:
            dupe_variant.delete()
            return 1
        except Exception as e:
            print(f"Could not delete {dupe_variant.pk}: {e}")
        return 0

    def _find_bad_symbolic_via_clinvar(self, dry_run: bool, genome_build, alt_seq):
        """ Due to VariantCoordinate.as_symbolic_variant bug - some historical data was incorrectly imported """

        fixed = 0

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

        for clinvar_variation_id, bad_variant in clinvar_variation_bad.items():
            original_variant = clinvar_variation_original[clinvar_variation_id]
            print(f"{bad_variant} should have been {original_variant}")
            if not dry_run:
                fixed += self._merge_variant_dupe(bad_variant, original_variant)

        print(f"{genome_build} deleted {fixed} {alt_seq}")
