from django.core.management.base import BaseCommand

from annotation.models import AnnotationRangeLock
from genes.hgvs import HGVSMatcher
from snpdb.models import Variant, Sequence, GenomeBuild, Locus
from django.db.models.functions import Length
from django.db.models import Q


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
        dupes = []
        dry_run = options["dry_run"]

        for genome_build in GenomeBuild.builds_with_annotation():
            matcher = HGVSMatcher(genome_build)

            q_contig = Variant.get_contigs_q(genome_build)
            for v in long_variants.filter(q_contig):
                vc = v.coordinate.as_internal_symbolic(genome_build)

                old_hgvs = matcher.variant_coordinate_to_g_hgvs(v.coordinate)
                new_hgvs = matcher.variant_coordinate_to_g_hgvs(vc)

                if old_hgvs != new_hgvs:
                    raise ValueError(f"Variant: {v.pk} old: {old_hgvs} != {new_hgvs}")

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
                            if v.cohortgenotype_set.exists():
                                dupes.append((v, existing))
                                continue
                            v.clinvar_set.all().update(variant=existing)
                            AnnotationRangeLock.release_variant(v)

                            # Classification protects Variant FK, so won't delete if that exists
                            v.delete()
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

        print(f"Non-symbolic not converted: {len(not_symbolic)}")
        print(f"Dupes: {len(dupes)}")
