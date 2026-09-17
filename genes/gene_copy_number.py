"""
Whole-gene copy number calls - "EGFR amplification" - as GeneCopyNumberEvent records.

@see snpdb.gene_level_variants for why one of these is a Variant, and genes.gene_level_resolver for
turning the caller's gene name into the identity it is stored under.

Identity is the gene plus the direction. The segment coordinates a caller writes are the panel's
target window rather than the event - they move with the manifest, differ between assays, and are
absent altogether from the classifications labs already hold - so they are per observation and never
part of what makes two calls the same thing.

'amplification' and 'loss' are the words written out; 'amp', 'gain', 'deletion' and 'del' are
accepted from a lab or a search box. 'deletion' stays input-only because as output it reads as a
coordinate event.
"""
import re
from dataclasses import dataclass
from typing import Optional

from genes.gene_level_resolver import GeneLevelNameResolver
from genes.models import (
    GeneCopyNumberEvent,
    GeneCopyNumberEventKind,
    GeneLevelId,
    gene_copy_number_canonical_str,
)
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from snpdb.gene_level_variants import (
    GENE_LEVEL_CONTIG_NAME,
    GENE_LEVEL_REF,
    GENE_LEVEL_SVLEN,
)
from snpdb.models import Variant, VariantCoordinate

# What a lab or a user writes after the gene, lower-cased
COPY_NUMBER_KIND_WORDS = {
    "amplification": GeneCopyNumberEventKind.GAIN,
    "amp": GeneCopyNumberEventKind.GAIN,
    "gain": GeneCopyNumberEventKind.GAIN,
    "loss": GeneCopyNumberEventKind.LOSS,
    "deletion": GeneCopyNumberEventKind.LOSS,
    "del": GeneCopyNumberEventKind.LOSS,
}
# 'EGFR amplification' - one gene then one word. Longest word first so 'amp' doesn't shadow
# 'amplification', and the space between them is optional because a classification's imported c.HGVS
# reaches us with its spaces removed (@see ImportedAlleleInfo._tidy_input_value). A gene that has to
# resolve is what keeps that from claiming ordinary strings - 'CAMP' asks for a gene named 'C'.
COPY_NUMBER_STRING_PATTERN = re.compile(
    rf"^\s*([A-Za-z0-9.\-]+?)\s*({'|'.join(sorted(COPY_NUMBER_KIND_WORDS, key=len, reverse=True))})\s*$",
    re.IGNORECASE)


@dataclass(frozen=True)
class ResolvedGeneCopyNumberEvent:
    """ A copy number event's identity, before it has a Variant. The gene is both the
        Locus.position and what the alt carries - @see snpdb.gene_level_variants """
    gene: GeneLevelId
    kind: GeneCopyNumberEventKind

    @property
    def alt(self) -> str:
        """ The alt repeats the gene the position already names, so the alt alone says what the
            variant is """
        alt_kind = GeneCopyNumberEventKind(self.kind).alt_kind
        return GeneLevelSymbolicAlt.format(alt_kind, self.gene.alt_namespace, self.gene.pk)

    @property
    def variant_coordinate(self) -> VariantCoordinate:
        """ What the loader writes as a VCF record, and looks the Variant up by afterwards """
        return VariantCoordinate(chrom=GENE_LEVEL_CONTIG_NAME, position=self.gene.pk,
                                 ref=GENE_LEVEL_REF, alt=self.alt, svlen=GENE_LEVEL_SVLEN)

    @property
    def canonical_str(self) -> str:
        """ The same string GeneCopyNumberEvent.canonical_str gives, before the Variant exists """
        return gene_copy_number_canonical_str(self.gene, self.kind)


def parse_gene_copy_number_string(copy_number_string: str) -> Optional[tuple[str, GeneCopyNumberEventKind]]:
    """ 'egfr amp' -> ('egfr', GAIN). The gene name is returned as written, for a resolver to place """
    if m := COPY_NUMBER_STRING_PATTERN.match(copy_number_string or ""):
        gene_name, word = m.groups()
        return gene_name, COPY_NUMBER_KIND_WORDS[word.lower()]
    return None


def resolve_gene_copy_number_string(copy_number_string: str,
                                    resolver: GeneLevelNameResolver = None) \
        -> Optional[ResolvedGeneCopyNumberEvent]:
    """ 'EGFR amplification' -> the identity it will be stored under, whose variant_coordinate goes
        through the VCF insert pipeline like any other coordinate.

        The gene has to be one we already know, so an arbitrary word pair doesn't mint an identity.
        @see ImportedAlleleInfo for where a classification target comes in this way. """

    parsed = parse_gene_copy_number_string(copy_number_string)
    if parsed is None:
        return None

    gene_name, kind = parsed
    if resolver is None:
        resolver = GeneLevelNameResolver()
    if resolved_gene := resolver.resolve_gene(gene_name, allow_unknown=False):
        return ResolvedGeneCopyNumberEvent(gene=resolved_gene.gene_level_id, kind=kind)
    return None


def find_gene_copy_number_events_for_string(copy_number_string: str,
                                            resolver: GeneLevelNameResolver = None) \
        -> list[GeneCopyNumberEvent]:
    """ Lookup only - 'EGFR amplification' finds the event if we have it, and mints nothing if we
        don't. Search runs on whatever a user types, so it must not create identities. """

    parsed = parse_gene_copy_number_string(copy_number_string)
    if parsed is None:
        return []

    gene_name, kind = parsed
    if resolver is None:
        resolver = GeneLevelNameResolver()
    symbols = {resolver.canonical_symbol(n) for n in resolver.split_gene_names(gene_name)}
    return list(GeneCopyNumberEvent.objects.filter(gene__symbol_str__in=symbols, kind=kind)
                .select_related("variant", "gene"))


def create_gene_copy_number_events_for_variants(variant_qs) -> int:
    """ The GeneCopyNumberEvent rows for gene-level variants an insert pipeline has just created.

        Everything an event holds is already in the Variant - the gene is Locus.position and the alt
        carries the direction - so this reads the variants rather than the file they came from, and
        one implementation serves every loader.

        :return: how many were created """

    events = []
    for variant in variant_qs.filter(Variant.get_gene_level_q(), genecopynumberevent__isnull=True) \
                             .select_related("locus", "alt"):
        alt_kind, _namespace, _gene_id = GeneLevelSymbolicAlt.parse(variant.alt.seq)
        kind = GeneCopyNumberEventKind.from_alt_kind(alt_kind)
        if kind is None:
            continue
        events.append(GeneCopyNumberEvent(variant=variant, gene_id=variant.locus.position, kind=kind))
    created = GeneCopyNumberEvent.objects.bulk_create(events, ignore_conflicts=True)
    return len(created)
