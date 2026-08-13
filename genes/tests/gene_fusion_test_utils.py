"""
Building a fusion Variant directly, for tests.

Production has one way in - the VCF insert pipeline (@see snpdb.gene_level_variants) - which is far
more machinery than a test wanting one fusion to exist needs. This makes the same rows that pipeline
would, from the same ResolvedFusion, so a test fixture cannot drift from what a real import produces.
"""
from typing import Optional

from django.db import transaction

from genes.gene_fusions import GeneFusionResolver, create_gene_fusions_for_variants
from genes.models import GeneFusion
from library.utils import sha256sum_str
from snpdb.gene_level_variants import GENE_LEVEL_REF, GENE_LEVEL_SVLEN
from snpdb.models import Contig, Locus, Sequence, Variant


def _get_sequence(seq: str) -> Sequence:
    sequence, _ = Sequence.objects.get_or_create(seq=seq, defaults={"seq_sha256_hash": sha256sum_str(seq)})
    return sequence


@transaction.atomic
def create_gene_fusion(gene_a: Optional[str], gene_b: Optional[str], directionality_known: bool = True,
                       resolver: GeneFusionResolver = None) -> GeneFusion:
    """ ('BCR', 'ABL1') -> the GeneFusion, as the insert pipeline would have created it """

    if resolver is None:
        resolver = GeneFusionResolver()
    resolved_fusion = resolver.resolve_fusion(resolver.resolve_side(gene_a) if gene_a else None,
                                              resolver.resolve_side(gene_b) if gene_b else None,
                                              directionality_known)
    variant_coordinate = resolved_fusion.variant_coordinate
    locus, _ = Locus.objects.get_or_create(contig=Contig.get_gene_level(),
                                           position=variant_coordinate.position,
                                           ref=_get_sequence(GENE_LEVEL_REF))
    variant, _ = Variant.objects.get_or_create(locus=locus, alt=_get_sequence(variant_coordinate.alt),
                                               svlen=GENE_LEVEL_SVLEN,
                                               defaults={"end": variant_coordinate.position})
    create_gene_fusions_for_variants(Variant.objects.filter(pk=variant.pk))
    return GeneFusion.objects.get(variant=variant)
