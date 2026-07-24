""" Q object builders for the standard "which variants?" filters - contigs, gene symbols and variant types.

    Shared by the All Variants page grid (variantopedia.grids.AllVariantsGrid) and the analysis pipeline's
    AllVariantsNode, so both compose the same queries.
"""
import operator
from functools import reduce
from typing import Any, Iterable, Optional

from django.conf import settings
from django.contrib.auth.models import User
from django.db.models import Q

from annotation.models import AnnotationVersion, VariantTranscriptAnnotation
from genes.models import Gene, GeneSymbol
from snpdb.models.models_enums import SequenceRole
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_user_settings import AllVariantsFilter
from snpdb.models.models_variant import Variant

# The smallest standard autosome - a brand new user's default, so their first page load is cheap
DEFAULT_CONTIG_NAME = "21"


class VariantType:
    """ Tokens for the variant type filter, stored in AllVariantsFilter.filters.
        Symbolic alts (eg '<DEL>') are used as tokens directly. """
    REFERENCE = "reference"
    SNV = "snv"
    INDEL = "indel"
    COMPLEX = "complex"
    SYMBOLIC = "symbolic"  # Any variant with an SVLEN


_VARIANT_TYPE_Q_FUNCS = {
    VariantType.REFERENCE: Variant.get_reference_q,
    VariantType.SNV: Variant.get_snp_q,
    VariantType.INDEL: Variant.get_indel_q,
    VariantType.COMPLEX: Variant.get_complex_subsitution_q,
    VariantType.SYMBOLIC: Variant.get_symbolic_q,
}

VARIANT_TYPE_LABELS = {
    VariantType.REFERENCE: "Reference",
    VariantType.SNV: "SNV",
    VariantType.INDEL: "Indel",
    VariantType.COMPLEX: "Complex sub",
    VariantType.SYMBOLIC: "Structural",
}

# The types offered on the All Variants page - reference variants are always excluded there, and the
# symbolic alts are broken out into a button each (appended by get_all_variant_types)
STANDARD_VARIANT_TYPES = [VariantType.SNV, VariantType.INDEL, VariantType.COMPLEX]


def get_symbolic_variant_types() -> list[str]:
    """ Symbolic alts (eg '<DEL>') this deployment accepts - empty when symbolic alts are disabled """
    if settings.VARIANT_SYMBOLIC_ALT_ENABLED:
        return sorted(settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES)
    return []


def get_all_variant_types() -> list[str]:
    return STANDARD_VARIANT_TYPES + get_symbolic_variant_types()


def get_variant_type_label(variant_type: str) -> str:
    """ Symbolic alts display as eg 'DEL' - the angle brackets are noise on a button """
    return VARIANT_TYPE_LABELS.get(variant_type) or variant_type.strip("<>")


def get_variant_type_q(variant_type: str) -> Q:
    if q_func := _VARIANT_TYPE_Q_FUNCS.get(variant_type):
        return q_func()
    if variant_type in settings.VARIANT_SYMBOLIC_ALT_VALID_TYPES:
        # Symbolic alts are only meaningful with an SVLEN - @see issue #1663
        return Q(alt__seq=variant_type) & Variant.get_symbolic_q()
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


def get_non_standard_contig_ids(genome_build: GenomeBuild) -> list[int]:
    """ Alt scaffolds, patches and unplaced/unlocalized contigs """
    qs = genome_build.contigs.exclude(role=SequenceRole.ASSEMBLED_MOLECULE)
    return list(qs.values_list("pk", flat=True))


def get_contigs_q(genome_build: GenomeBuild, contig_ids: Optional[Iterable[int]] = None,
                  non_standard_contigs: bool = False) -> Q:
    """ Restrict to the build's contigs, narrowed to a contig selection when there is one """
    contig_id_list = list(contig_ids or [])
    if non_standard_contigs:
        contig_id_list.extend(get_non_standard_contig_ids(genome_build))
    if contig_id_list:
        return Q(locus__contig_id__in=contig_id_list)
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


def get_gene_symbols_q(annotation_version: AnnotationVersion, gene_symbols: Optional[Iterable[Any]],
                       traverse_aliases: bool = True) -> Optional[Q]:
    symbols = resolve_gene_symbols(gene_symbols)
    if not symbols:
        return None
    genes = get_genes_for_gene_symbols(symbols, traverse_aliases=traverse_aliases)
    # pk__in form so a variant overlapping several of the genes still returns a single row
    return VariantTranscriptAnnotation.get_overlapping_genes_q(annotation_version.variant_annotation_version, genes)


def get_variant_filter_q(genome_build: GenomeBuild, annotation_version: AnnotationVersion, *,
                         contig_ids: Optional[Iterable[int]] = None, non_standard_contigs: bool = False,
                         gene_symbols: Optional[Iterable[Any]] = None,
                         variant_types: Optional[Iterable[str]] = None) -> Q:
    """ The standard variant filters composed into a single Q """
    filter_list = [get_contigs_q(genome_build, contig_ids=contig_ids, non_standard_contigs=non_standard_contigs)]
    if (q_genes := get_gene_symbols_q(annotation_version, gene_symbols)) is not None:
        filter_list.append(q_genes)
    if (q_types := get_variant_types_q(variant_types)) is not None:
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
