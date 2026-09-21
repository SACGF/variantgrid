""" Q object builders for the standard "which variants?" filters - contigs, gene symbols and variant types.

    Shared by the All Variants page grid (variantopedia.grids.AllVariantsGrid) and the analysis pipeline's
    AllVariantsNode, so both compose the same queries.
"""
import operator
from collections.abc import Iterable
from functools import lru_cache, reduce
from typing import Any, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Max, Min, Q

from annotation.models import AnnotationVersion, VariantGeneOverlap, VariantTranscriptAnnotation
from genes.models import Gene, GeneSymbol
from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from snpdb.models.models_enums import SequenceRole
from snpdb.models.models_genome import Contig, GenomeBuild
from snpdb.models.models_user_settings import AllVariantsFilter
from snpdb.models.models_variant import Sequence, Variant

# The smallest standard autosome - a brand new user's default, so their first page load is cheap
DEFAULT_CONTIG_NAME = "21"


class VariantType:
    """ Tokens for the variant type filter, stored in AllVariantsFilter.filters.
        Symbolic alts (eg '<DEL>') are used as tokens directly. """
    REFERENCE = "reference"
    SNV = "snv"
    INDEL = "indel"
    COMPLEX = "complex"
    SYMBOLIC = "symbolic"  # Structural variants - a symbolic alt with coordinates, ie not gene-level
    # Gene-level events - @see snpdb.gene_level_variants
    FUSION = "fusion"
    COPY_NUMBER = "copy_number"
    SPLICE = "splice"


GENE_LEVEL_VARIANT_TYPES = [VariantType.FUSION, VariantType.COPY_NUMBER, VariantType.SPLICE]


def _lookup_sequence_ids() -> dict[str, int]:
    return Sequence.get_pk_by_seq(Q(seq__in=[*Variant.BASES, Variant.REFERENCE_ALT]))


@lru_cache
def _cached_sequence_ids() -> dict[str, int]:
    return _lookup_sequence_ids()


def _get_sequence_ids() -> dict[str, int]:
    """ pks of the single base sequences and of the reference alt - the type filters compare these
        against Variant.alt_id / Locus.ref_id rather than joining snpdb_sequence to test seq, which
        collapses the planner's row estimate and costs the streaming plan (#1887).

        Sequence rows are never deleted, so the pks are cached - except under test, where each
        database builds its own (@see library.guardian_utils.admin_bot). A database that doesn't have
        them all yet (nothing imported) is looked up again rather than cached as missing """
    if settings.UNIT_TEST:
        return _lookup_sequence_ids()
    sequence_ids = _cached_sequence_ids()
    if len(sequence_ids) <= len(Variant.BASES):
        _cached_sequence_ids.cache_clear()
        sequence_ids = _lookup_sequence_ids()
    return sequence_ids


def _base_sequence_ids() -> list[int]:
    sequence_ids = _get_sequence_ids()
    return [pk for seq, pk in sequence_ids.items() if seq in Variant.BASES]


def _base_and_reference_sequence_ids() -> list[int]:
    return list(_get_sequence_ids().values())


def get_snv_q() -> Q:
    base_ids = _base_sequence_ids()
    return Q(locus__ref_id__in=base_ids, alt_id__in=base_ids)


def _get_plain_q() -> Q:
    """ Not symbolic and not gene-level - both of those carry an SVLEN """
    return Q(svlen__isnull=True)


def get_indel_q() -> Q:
    """ One side a single base, the other a longer sequence """
    base_ids = _base_sequence_ids()
    other_ids = _base_and_reference_sequence_ids()
    insertion = Q(locus__ref_id__in=base_ids) & ~Q(alt_id__in=other_ids)
    deletion = Q(alt_id__in=base_ids) & ~Q(locus__ref_id__in=other_ids)
    return _get_plain_q() & (insertion | deletion)


def get_complex_substitution_q() -> Q:
    """ Both sides longer than a single base """
    other_ids = _base_and_reference_sequence_ids()
    return _get_plain_q() & ~Q(locus__ref_id__in=other_ids) & ~Q(alt_id__in=other_ids)


def _alt_in_q(q_sequence: Q) -> Q:
    """ An alt_id IN (subquery) rather than a join to Sequence - OR'd with the other types, a join
        stops Postgres using the per-branch plans and the whole-build scan goes from ~1s to ~14s """
    return Q(alt__in=Sequence.objects.filter(q_sequence))


def _gene_level_kinds_q(*kinds: str) -> Q:
    """ By alt prefix - '<FUSION' covers FUSION_UNORDERED too """
    q_sequence = reduce(operator.or_, [Q(seq__startswith=f"<{kind}") for kind in kinds])
    return Variant.get_gene_level_q() & _alt_in_q(q_sequence)


def get_structural_variant_q() -> Q:
    """ Symbolic alts with coordinates - gene-level events also have an SVLEN (0), so name the alts """
    return Variant.get_symbolic_q() & _alt_in_q(Q(seq__in=settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES))


_VARIANT_TYPE_Q_FUNCS = {
    VariantType.REFERENCE: Variant.get_reference_q,
    VariantType.SNV: get_snv_q,
    VariantType.INDEL: get_indel_q,
    VariantType.COMPLEX: get_complex_substitution_q,
    VariantType.SYMBOLIC: get_structural_variant_q,
    VariantType.FUSION: lambda: _gene_level_kinds_q(GeneLevelSymbolicAlt.FUSION),
    VariantType.COPY_NUMBER: lambda: _gene_level_kinds_q(GeneLevelSymbolicAlt.GAIN, GeneLevelSymbolicAlt.LOSS),
    VariantType.SPLICE: lambda: _gene_level_kinds_q(GeneLevelSymbolicAlt.SPLICE),
}

VARIANT_TYPE_LABELS = {
    VariantType.REFERENCE: "Reference",
    VariantType.SNV: "SNV",
    VariantType.INDEL: "Indel",
    VariantType.COMPLEX: "Complex sub",
    VariantType.SYMBOLIC: "Structural",
    VariantType.FUSION: "Fusion",
    VariantType.COPY_NUMBER: "Copy number",
    VariantType.SPLICE: "Splicing",
}

# The types offered on the All Variants page - reference variants are always excluded there, and the
# symbolic alts are broken out into a button each (appended by get_all_variant_types)
STANDARD_VARIANT_TYPES = [VariantType.SNV, VariantType.INDEL, VariantType.COMPLEX]


def get_symbolic_variant_types() -> list[str]:
    """ Symbolic alts (eg '<DEL>') this deployment accepts - empty when symbolic alts are disabled """
    if settings.VARIANT_SYMBOLIC_ALT_ENABLED:
        return sorted(settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES)
    return []


def get_gene_level_variant_types() -> list[str]:
    """ Empty when gene-level variants are disabled """
    if settings.VARIANT_GENE_LEVEL_ENABLED:
        return list(GENE_LEVEL_VARIANT_TYPES)
    return []


def get_all_variant_types() -> list[str]:
    return STANDARD_VARIANT_TYPES + get_symbolic_variant_types() + get_gene_level_variant_types()


def get_variant_type_label(variant_type: str) -> str:
    """ Symbolic alts display as eg 'DEL' - the angle brackets are noise on a button """
    return VARIANT_TYPE_LABELS.get(variant_type) or variant_type.strip("<>")


def get_variant_type_q(variant_type: str) -> Q:
    if q_func := _VARIANT_TYPE_Q_FUNCS.get(variant_type):
        return q_func()
    if variant_type in settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES:
        # Symbolic alts are only meaningful with an SVLEN - @see issue #1663
        return _alt_in_q(Q(seq=variant_type)) & Variant.get_symbolic_q()
    msg = f"Unknown variant type: '{variant_type}'"
    raise ValueError(msg)


def get_variant_types_q(variant_types: Optional[Iterable[str]],
                        all_variant_types: Optional[Iterable[str]] = None) -> Optional[Q]:
    """ Returns None (ie no restriction) when the selection is empty or covers every available type """
    if not variant_types:
        return None
    selected = set(variant_types)
    if all_variant_types is None:
        all_variant_types = get_all_variant_types()
    if selected.issuperset(all_variant_types):
        return None
    return reduce(operator.or_, [get_variant_type_q(vt) for vt in sorted(selected)])


def get_gene_level_contig_ids(genome_build: GenomeBuild) -> list[int]:
    """ The shared coordinate-free contig fusions live on - @see snpdb.gene_level_variants """
    qs = genome_build.contigs.filter(role=SequenceRole.VG_GENE_LEVEL_FAKE_CONTIG)
    return list(qs.values_list("pk", flat=True))


def get_non_standard_contig_ids(genome_build: GenomeBuild) -> list[int]:
    """ Alt scaffolds, patches and unplaced/unlocalized contigs """
    qs = genome_build.contigs.exclude(role=SequenceRole.ASSEMBLED_MOLECULE)
    return list(qs.values_list("pk", flat=True))


def get_contig_ids_for_variant_types(genome_build: GenomeBuild, contig_ids: Iterable[int],
                                     variant_types: Optional[Iterable[str]]) -> list[int]:
    """ The gene-level contig is not a chromosome anyone can tick, so a gene-level type is what lets
        it through - without it a contig selection hides every fusion whatever else is on """
    contig_id_list = list(contig_ids)
    if contig_id_list and (variant_types is None or set(variant_types) & set(GENE_LEVEL_VARIANT_TYPES)):
        contig_id_list.extend(get_gene_level_contig_ids(genome_build))
    return contig_id_list


def get_contigs_q(genome_build: GenomeBuild, contig_ids: Optional[Iterable[int]] = None,
                  non_standard_contigs: bool = False) -> Q:
    """ Restrict to the build's contigs, narrowed to a contig selection when there is one """
    contig_id_list = list(contig_ids or [])
    if non_standard_contigs:
        contig_id_list.extend(get_non_standard_contig_ids(genome_build))
    if contig_id_list:
        # A gene's contigs can arrive twice (its own and the gene-level one) - keep the IN list unique
        return Q(locus__contig_id__in=sorted(set(contig_id_list)))
    return Variant.get_contigs_q(genome_build)


def resolve_gene_symbols(gene_symbols: Optional[Iterable[Any]]) -> list[GeneSymbol]:
    """ Accepts GeneSymbol instances or symbol strings (GeneSymbol's pk is the symbol) """
    symbols: list[GeneSymbol] = []
    symbol_strs: list[str] = []
    for gene_symbol in gene_symbols or []:
        if isinstance(gene_symbol, GeneSymbol):
            symbols.append(gene_symbol)
        else:
            symbol_strs.append(gene_symbol)
    if symbol_strs:
        symbols.extend(GeneSymbol.objects.filter(pk__in=symbol_strs))
    return symbols


def get_genes_for_gene_symbols(gene_symbols: Iterable[GeneSymbol], traverse_aliases: bool = True) -> set[Gene]:
    genes: set[Gene] = set()
    for gene_symbol in gene_symbols:
        if traverse_aliases:
            genes |= gene_symbol.alias_meta.genes
        else:
            genes |= set(gene_symbol.genes)
    return genes


def get_gene_symbol_alias_strs(gene_symbol: GeneSymbol) -> list[str]:
    """ The symbol plus every alias that resolves to the same genes """
    return gene_symbol.alias_meta.alias_symbol_strs


def get_gene_bounds_q(annotation_version: AnnotationVersion, genes: Iterable[Gene]) -> Optional[Q]:
    """ The locus range the genes' variants span, per gene and contig.

        The overlap Q alone is a pk IN (subquery), so a gene filter walks the build in genomic order
        testing every locus until it reaches the gene. AND-ing the bounds on lets the locus index seek
        straight there - they come off the same rows the overlap Q matches, so nothing new is excluded.
        None means the genes overlap no variants in this annotation version (#1887) """
    qs = VariantGeneOverlap.objects.filter(version=annotation_version.variant_annotation_version, gene__in=genes)
    bounds = qs.values("gene_id", "variant__locus__contig_id") \
        .annotate(min_position=Min("variant__locus__position"), max_position=Max("variant__locus__position"))
    q_list = [Q(locus__contig_id=b["variant__locus__contig_id"],
                locus__position__range=(b["min_position"], b["max_position"])) for b in bounds]
    if not q_list:
        return None
    return reduce(operator.or_, q_list)


def get_gene_symbols_q(annotation_version: AnnotationVersion, gene_symbols: Optional[Iterable[Any]],
                       traverse_aliases: bool = True) -> Optional[Q]:
    symbols = resolve_gene_symbols(gene_symbols)
    if not symbols:
        return None
    genes = get_genes_for_gene_symbols(symbols, traverse_aliases=traverse_aliases)
    # pk__in form so a variant overlapping several of the genes still returns a single row
    q_overlap = VariantTranscriptAnnotation.get_overlapping_genes_q(
        annotation_version.variant_annotation_version, genes)
    q_bounds = get_gene_bounds_q(annotation_version, genes)
    if q_bounds is None:
        return Q(pk__isnull=True)  # No variants overlap the genes - the answer is already known
    return q_overlap & q_bounds


def get_contig_ids_for_gene_symbols(genome_build: GenomeBuild, gene_symbols: Optional[Iterable[Any]],
                                    traverse_aliases: bool = True) -> list[int]:
    """ The contigs a gene's transcripts live on - a gene filter is only visible once its contig is selected """
    symbols = resolve_gene_symbols(gene_symbols)
    if not symbols:
        return []
    genes = get_genes_for_gene_symbols(symbols, traverse_aliases=traverse_aliases)
    # TranscriptVersion.contig exists as an optimisation to restrict Variant queries
    contig_qs = Contig.objects.filter(transcriptversion__genome_build=genome_build,
                                      transcriptversion__gene_version__gene__in=genes)
    return sorted(set(contig_qs.values_list("pk", flat=True)))


def get_variant_filter_q(genome_build: GenomeBuild, annotation_version: AnnotationVersion, *,
                         contig_ids: Optional[Iterable[int]] = None, non_standard_contigs: bool = False,
                         gene_symbols: Optional[Iterable[Any]] = None,
                         variant_types: Optional[Iterable[str]] = None) -> Q:
    """ The standard variant filters composed into a single Q """
    contig_id_list = list(contig_ids or [])
    variant_type_list = list(variant_types) if variant_types is not None else None
    if gene_symbols:
        # A gene selection is a contig restriction of its own - the genes' contigs (so a gene on a
        # chromosome the user hasn't ticked still shows) plus the gene-level contig its fusions are on.
        # Without it a gene with no chromosome chosen scans the whole build instead of the gene bounds.
        contig_id_list.extend(get_contig_ids_for_gene_symbols(genome_build, gene_symbols))
        contig_id_list.extend(get_gene_level_contig_ids(genome_build))
    if contig_id_list or non_standard_contigs:
        contig_id_list = get_contig_ids_for_variant_types(genome_build, contig_id_list, variant_type_list)

    filter_list = [get_contigs_q(genome_build, contig_ids=contig_id_list,
                                 non_standard_contigs=non_standard_contigs)]
    if (q_genes := get_gene_symbols_q(annotation_version, gene_symbols)) is not None:
        filter_list.append(q_genes)
    if (q_types := get_variant_types_q(variant_type_list)) is not None:
        filter_list.append(q_types)
    return reduce(operator.and_, filter_list)


def get_default_all_variants_filters(genome_build: GenomeBuild) -> dict:
    """ What a user sees on the All Variants page before they've chosen anything """
    contig_ids = list(genome_build.standard_contigs.filter(name=DEFAULT_CONTIG_NAME).values_list("pk", flat=True))
    return {
        "contig_ids": contig_ids,
        "non_standard_contigs": False,
        "gene_symbols": [],
        "variant_types": get_all_variant_types(),
        "min_count": 0,
    }


def get_all_variants_filters(user: User, genome_build: GenomeBuild) -> dict:
    """ A user's saved All Variants page filters, falling back to the defaults """
    all_variants_filter = AllVariantsFilter.get(user, genome_build)
    return all_variants_filter.filters or get_default_all_variants_filters(genome_build)


def is_selective(filters: dict) -> bool:
    """ The All Variants page requires a filter that meaningfully restricts the scan - variant type and
        min count alone still walk the whole variant table """
    return bool(filters.get("contig_ids") or filters.get("non_standard_contigs") or filters.get("gene_symbols"))
