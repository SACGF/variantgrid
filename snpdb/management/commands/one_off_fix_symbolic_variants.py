from django.core.management.base import BaseCommand

from annotation.models import AnnotationRangeLock
from snpdb.models import Variant, Sequence, GenomeBuild, Locus
from django.db.models.functions import Length
from django.db.models import Q


class Command(BaseCommand):
    """
        There was a bug where SVs were not converted to symbolic for a while

        This will try and fix them
    """
    def handle(self, *args, **options):
        long_sequences = Sequence.objects.all().annotate(seq_length=Length("seq")).filter(seq_length__gte=1000)
        Variant.objects.filter(locus__ref=Q())
        long_variants = Variant.objects.filter(Q(locus__ref__in=long_sequences) | Q(alt__in=long_sequences))

        base_lookup = {s: Sequence.objects.get(seq=s) for s in ["G", "A", "T", "C", "<DEL>", "<DUP>"]}
        not_symbolic = []
        dupes = []

        for genome_build in GenomeBuild.builds_with_annotation():
            q_contig = Variant.get_contigs_q(genome_build)
            for v in long_variants.filter(q_contig):
                vc = v.coordinate.as_internal_symbolic()
                if not vc.is_symbolic():
                    not_symbolic.append(v)
                    continue

                try:
                    existing = Variant.get_from_variant_coordinate(vc, genome_build)

                    # I think these came in via ClinVar import - so will try to get rid of it, if it isn't referenced

                    # by anything else
                    if v.cohortgenotype_set.exists():
                        dupes.append((v, existing))
                        continue
                    v.clinvar_set.all().update(variant=existing)
                    AnnotationRangeLock.release_variant(v)

                    # Classification protects Variant FK, so won't delete if that exists
                    v.delete()
                except Variant.DoesNotExist:
                    print(f"Fixing {v}")
                    new_ref = base_lookup[vc.ref]
                    v.locus = Locus.objects.get_or_create(contig=v.locus.contig, position=vc.position, ref=new_ref)[0]
                    v.alt = base_lookup[vc.alt]
                    v.svlen = vc.svlen
                    v.save()

        print(f"Non-symbolic not converted: {len(not_symbolic)}")
        print(f"Dupes: {len(dupes)}")
