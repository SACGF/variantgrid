"""
Building a gene-level Variant directly, for tests.

Production has one way in - the VCF insert pipeline (@see snpdb.gene_level_variants) - which is far
more machinery than a test wanting one event to exist needs. This makes the same rows that pipeline
would, from the same resolved identity, so a test fixture cannot drift from what a real import
produces. Fusions build on this in gene_fusion_test_utils.
"""
from django.db import transaction

from genes.gene_copy_number import (
    ResolvedGeneCopyNumberEvent,
    create_gene_copy_number_events_for_variants,
)
from genes.gene_level_resolver import GeneLevelNameResolver
from genes.models import GeneCopyNumberEvent, GeneCopyNumberEventKind
from library.utils import sha256sum_str
from snpdb.gene_level_variants import GENE_LEVEL_SVLEN
from snpdb.models import Contig, Locus, Sequence, Variant, VariantCoordinate


def get_sequence(seq: str) -> Sequence:
    sequence, _ = Sequence.objects.get_or_create(seq=seq, defaults={"seq_sha256_hash": sha256sum_str(seq)})
    return sequence


def create_gene_level_variant(variant_coordinate: VariantCoordinate) -> Variant:
    """ The Variant the insert pipeline would have created for a resolved gene-level identity """
    locus, _ = Locus.objects.get_or_create(contig=Contig.get_gene_level(),
                                           position=variant_coordinate.position,
                                           ref=get_sequence(variant_coordinate.ref))
    variant, _ = Variant.objects.get_or_create(locus=locus, alt=get_sequence(variant_coordinate.alt),
                                               svlen=GENE_LEVEL_SVLEN,
                                               defaults={"end": variant_coordinate.position})
    return variant


@transaction.atomic
def create_gene_copy_number_event(gene_name: str, kind: GeneCopyNumberEventKind,
                                  resolver: GeneLevelNameResolver = None) -> GeneCopyNumberEvent:
    """ ('EGFR', GAIN) -> the GeneCopyNumberEvent, as the insert pipeline would have created it """

    if resolver is None:
        resolver = GeneLevelNameResolver()
    resolved_gene = resolver.resolve_gene(gene_name)
    event = ResolvedGeneCopyNumberEvent(gene=resolved_gene.gene_level_id, kind=kind)
    variant = create_gene_level_variant(event.variant_coordinate)
    create_gene_copy_number_events_for_variants(Variant.objects.filter(pk=variant.pk))
    return GeneCopyNumberEvent.objects.get(variant=variant)
