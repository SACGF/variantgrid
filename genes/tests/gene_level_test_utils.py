"""
Building a gene-level Variant directly, for tests.

Production has one way in - the VCF insert pipeline (@see snpdb.gene_level_variants) - which is far
more machinery than a test wanting one event to exist needs. This makes the same rows that pipeline
would, from the same resolved identity, so a test fixture cannot drift from what a real import
produces. Fusions build on this in gene_fusion_test_utils.
"""
from django.db import transaction

from annotation.tests.test_data_fake_genes import _create_fake_gene_version, _insert_transcript_data
from genes.gene_copy_number import (
    ResolvedGeneCopyNumberEvent,
    create_gene_copy_number_events_for_variants,
)
from genes.gene_level_resolver import GeneLevelNameResolver
from genes.gene_splice import (
    ResolvedSpliceEvent,
    SpliceEventVariant,
    canonical_splice_label,
)
from genes.models import GeneCopyNumberEvent, GeneCopyNumberEventKind, ReleaseTranscriptVersion
from genes.models_enums import AnnotationConsortium
from library.utils import sha256sum_str
from snpdb.gene_level_variants import GENE_LEVEL_SVLEN
from snpdb.models import Contig, Locus, Sequence, Variant, VariantCoordinate


def get_sequence(seq: str) -> Sequence:
    """ Upper-cased, as every path that inserts a Sequence in production does """
    seq = seq.upper()
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


@transaction.atomic
def create_splice_event_variant(gene_name: str, label: str,
                                resolver: GeneLevelNameResolver = None) -> SpliceEventVariant:
    """ ('AR', 'V7') -> the splice Variant, as the insert pipeline would have created it. The label
        is canonicalised on the way in, as an import's is, so a fixture holds what a real one does.
        A splice event has no row of its own, so this hands back the view of the Variant """

    canonical = canonical_splice_label(label)
    if canonical is None:
        raise ValueError(f"'{label}' is not a splice label shape we accept")
    if resolver is None:
        resolver = GeneLevelNameResolver()
    resolved_gene = resolver.resolve_gene(gene_name)
    event = ResolvedSpliceEvent(gene=resolved_gene.gene_level_id, label=canonical)
    variant = create_gene_level_variant(event.variant_coordinate)
    return SpliceEventVariant(variant=variant, gene=event.gene, label=canonical)


def make_release_gene(genome_build, release, gene_id, gene_symbol, transcript_id, contig, start,
                      hgnc_id=None):
    """ A gene of the release with one 10 kb transcript from start, so a position inside it resolves """
    gene_version = _create_fake_gene_version(genome_build, gene_id, gene_symbol,
                                             AnnotationConsortium.ENSEMBL)
    gene_version.hgnc_id = hgnc_id
    gene_version.save()
    data = {
        "id": transcript_id,
        "gene_name": gene_symbol,
        "biotype": [],
        "genome_builds": {
            genome_build.name: {
                "url": "fake",
                "exons": [[start, start + 10_000, 0, 1, 10_001, None]],
                "contig": contig,
                "strand": "+",
                "cds_end": start + 10_000,
                "cds_start": start,
            }
        },
    }
    transcript_version = _insert_transcript_data(genome_build, data, gene_version, release)
    ReleaseTranscriptVersion.objects.get_or_create(release=release, transcript_version=transcript_version)
    return gene_version.gene
