"""
Building a fusion Variant directly, for tests.

Production has one way in - the VCF insert pipeline (@see snpdb.gene_level_variants) - which is far
more machinery than a test wanting one fusion to exist needs. This makes the same rows that pipeline
would, from the same ResolvedFusion, so a test fixture cannot drift from what a real import produces.
"""
from typing import Optional

from django.db import transaction

from genes.gene_fusions import GeneFusionResolver, ResolvedFusion, create_gene_fusions_for_variants
from genes.models import GeneFusion, GeneLevelId
from genes.tests.gene_level_test_utils import create_gene_level_variant
from snpdb.models import Variant


@transaction.atomic
def create_gene_fusion(gene_a: Optional[str], gene_b: Optional[str], directionality_known: bool = True,
                       resolver: GeneFusionResolver = None) -> GeneFusion:
    """ ('BCR', 'ABL1') -> the GeneFusion, as the insert pipeline would have created it """

    if resolver is None:
        resolver = GeneFusionResolver()
    resolved_fusion = resolver.resolve_fusion(resolver.resolve_side(gene_a) if gene_a else None,
                                              resolver.resolve_side(gene_b) if gene_b else None,
                                              directionality_known)
    return _create_from_resolved_fusion(resolved_fusion)


@transaction.atomic
def create_gene_fusion_for_ids(anchor: GeneLevelId, partner: Optional[GeneLevelId] = None,
                               is_ordered: bool = True) -> GeneFusion:
    """ The same rows, from identities a test built itself - for a side whose resolution is the
        thing under test rather than the thing being set up """
    return _create_from_resolved_fusion(ResolvedFusion(anchor=anchor, partner=partner,
                                                       is_ordered=is_ordered))


def _create_from_resolved_fusion(resolved_fusion: ResolvedFusion) -> GeneFusion:
    variant = create_gene_level_variant(resolved_fusion.variant_coordinate)
    create_gene_fusions_for_variants(Variant.objects.filter(pk=variant.pk))
    return GeneFusion.objects.get(variant=variant)
