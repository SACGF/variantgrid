"""
Matching text to genes: GeneSymbolMatcher resolves gene names to GeneSymbol rows a batch at a time with
GeneSymbol.resolve_names (the rules GeneSymbol.resolve_name applies to one name), ReleaseGeneMatcher
records which Gene a symbol meant in a particular GeneAnnotationRelease (symbols
move between genes over time, so matching is per release), and tokenize_gene_symbols /
partition_oversized_names split gene-list text and reject tokens over MAX_GENE_SYMBOL_LENGTH with a
recorded warning.
"""
import logging
import re
from collections import defaultdict, namedtuple
from collections.abc import Iterable
from functools import cached_property
from typing import Optional

from django.db.models import Q
from django.db.models.functions import Upper

from genes.models import (
    GeneAnnotationRelease,
    GeneListGeneSymbol,
    GeneSymbol,
    GeneSymbolAlias,
    GeneVersion,
    ReleaseGeneSymbol,
    ReleaseGeneSymbolGene,
)
from genes.models.models_gene import GeneNameResolution
from library.utils import clean_string

TokenizeResult = namedtuple("TokenizeResult", ["valid", "oversized"])


class GeneSymbolMatcher:
    """ Resolves gene names with GeneSymbol.resolve_names, remembering each one so a matcher can serve a whole
        import. Resolve a file's names together up front: get_gene_symbol_id on a name not seen yet queries for it """

    def __init__(self):
        self._resolutions: dict[Optional[str], GeneNameResolution] = {}

    def resolve_names(self, names: Iterable[Optional[str]]) -> dict[Optional[str], GeneNameResolution]:
        names = set(names)
        if unresolved := names - self._resolutions.keys():
            self._resolutions.update(GeneSymbol.resolve_names(unresolved))
        return {name: self._resolutions[name] for name in names}

    def get_gene_symbol_id(self, original_gene_symbol: str) -> Optional[str]:
        return self.resolve_names([original_gene_symbol])[original_gene_symbol].gene_symbol_id

    @cached_property
    def _release_gene_matchers(self):
        return [ReleaseGeneMatcher(release) for release in GeneAnnotationRelease.objects.all()]

    def _match_symbols_to_genes_in_releases(self):
        for gm in self._release_gene_matchers:
            gm.match_unmatched_in_hgnc_and_gene_lists()

    def _match_unmatched_symbols_in_releases(self, gene_symbol_ids):
        for gm in self._release_gene_matchers:
            gm.match_unmatched_symbols(gene_symbol_ids)

    def create_gene_list_gene_symbols_from_text(self, gene_list, gene_text, save=True):
        tokens = tokenize_gene_symbols(gene_text)
        return self.create_gene_list_gene_symbols(gene_list, tokens.valid, save=save,
                                                  oversized_skipped=tokens.oversized)

    def create_gene_list_gene_symbols(self, gene_list, gene_names_list, modification_info=None, save=True,
                                      oversized_skipped=None):
        """ save=False only returns unsaved objects.
            oversized_skipped: names already filtered out upstream (e.g. by tokenize_gene_symbols) that
            should be reported alongside any oversized names found here. """
        valid_names, oversized_here = partition_oversized_names(gene_names_list)
        all_oversized = list(oversized_skipped or []) + oversized_here

        resolutions = self.resolve_names(valid_names)
        gene_list_gene_symbols = []
        ambiguous_aliases = []
        for original_name in valid_names:
            resolution = resolutions[original_name]
            if resolution.ambiguous_alias:
                ambiguous_aliases.append(str(resolution.ambiguous_alias))
            gene_list_gene_symbols.append(GeneListGeneSymbol(gene_list=gene_list,
                                                             original_name=original_name,
                                                             gene_symbol_id=resolution.gene_symbol_id,
                                                             gene_symbol_alias_id=resolution.gene_symbol_alias_id,
                                                             modification_info=modification_info))

        if save:
            GeneListGeneSymbol.objects.bulk_create(gene_list_gene_symbols, ignore_conflicts=True)
            if all_oversized:
                append_gene_list_warning(gene_list, format_oversized_warning(all_oversized))
            for ambiguous_alias in sorted(set(ambiguous_aliases)):
                append_gene_list_warning(gene_list, ambiguous_alias)
            gene_list.set_modified_to_now()
            # Only the symbols this list added can be unmatched - rematching every HGNC and gene list
            # symbol per release scans the whole gene list table and is what the import commands do
            gene_symbol_ids = {glgs.gene_symbol_id for glgs in gene_list_gene_symbols if glgs.gene_symbol_id}
            self._match_unmatched_symbols_in_releases(gene_symbol_ids)

        return gene_list_gene_symbols


class ReleaseGeneMatcher:
    """ Genes and symbols change over time. This creates DB records connecting them for a particular
        GeneAnnotationRelease.

        To actually retrieve matched genes for a release, use GeneAnnotationRelease.genes_for_symbol(s) """
    def __init__(self, release: GeneAnnotationRelease):
        self.release = release

    @cached_property
    def genes(self) -> dict[str, list]:
        gv_qs = GeneVersion.objects.filter(releasegeneversion__release=self.release).annotate(symbol_upper=Upper("gene_symbol"))
        genes_dict = defaultdict(list)
        for symbol_upper, gene_id in gv_qs.values_list("symbol_upper", "gene__identifier"):
            genes_dict[symbol_upper].append(gene_id)
        return genes_dict

    def _get_genes_dict(self):
        """
            We want to get all the ways a symbol can match a gene not through the release version
        """
        release_gene_version_symbols = defaultdict(set)
        # All the gene versions corresponding to genes in release
        release_gv_qs = GeneVersion.objects.filter(releasegeneversion__release=self.release)
        for gene_id, gene_symbol_id in release_gv_qs.values_list("gene_id", "gene_symbol_id"):
            release_gene_version_symbols[gene_id].add(gene_symbol_id)

        gene_versions_qs = GeneVersion.objects.filter(gene__in=self.release.get_genes()).select_related("genome_build")
        gene_versions_qs = gene_versions_qs.filter(gene_symbol__isnull=False)
        values = gene_versions_qs.values_list("gene_symbol_id", "gene_id", "version", "genome_build__name")
        genes_dict = defaultdict(dict)
        empty_set = set()
        for gene_symbol_id, gene_id, version, genome_build_name in values:
            if gene_symbol_id not in release_gene_version_symbols.get(gene_id, empty_set):
                match_info = f"Gene v{version}/{genome_build_name}"
                genes_dict[gene_symbol_id.upper()][gene_id] = match_info
        return genes_dict

    @cached_property
    def aliases_dict(self) -> dict[str, dict]:
        """ Upper case symbol -> {gene_id: match_info} for symbols the release doesn't have, reached either
            from an older GeneVersion's symbol or by a single GeneSymbolAlias hop (in either direction) to a
            symbol the release does have.

            One hop is deliberate: HGNC lists every previous symbol of a gene, so renames never need
            chaining, while chaining lets an alias string shared by two unrelated genes bridge them - e.g.
            'RP8' is an alias of both MT-TS2 and PDCD2, which used to match MT-TS2 gene lists to PDCD2 (#1669) """
        genes_dict = self._get_genes_dict()

        # Ordered so the match_info a symbol/gene pair reached by several alias rows gets is stable across runs
        alias_values = GeneSymbolAlias.objects.order_by("pk").values_list("alias", "gene_symbol_id", "source")
        for alias, symbol, source in alias_values.iterator():
            uc_alias = alias.upper()
            uc_symbol = symbol.upper()
            if uc_alias == uc_symbol:
                continue

            # "alias is an alias for symbol" - a gene list can use either end and mean the other
            for query, target_symbol in ((uc_alias, uc_symbol), (uc_symbol, uc_alias)):
                if query in self.genes:
                    continue  # Release has the symbol itself, so aliases are never consulted for it
                if gene_id_list := self.genes.get(target_symbol):
                    match_info = GeneSymbolAlias(alias=alias, gene_symbol_id=symbol, source=source).match_info
                    for gene_id in gene_id_list:
                        genes_dict[query].setdefault(gene_id, match_info)

        return genes_dict

    def _get_gene_id_and_match_info_for_symbol(self, gene_symbols) -> dict[str, list]:
        gene_symbol_gene_id_and_match_info = defaultdict(list)  # list items = (gene_id, match_info)
        for gene_symbol_id in gene_symbols:
            gene_name = clean_string(str(gene_symbol_id)).upper()

            if gene_id_list := self.genes.get(gene_name):
                for gene_id in gene_id_list:
                    gene_symbol_gene_id_and_match_info[gene_symbol_id].append((gene_id, None))
            elif alias_items := self.aliases_dict.get(gene_name):
                for gene_id, match_info in alias_items.items():
                    gene_symbol_gene_id_and_match_info[gene_symbol_id].append((gene_id, match_info))
            # Else - no match?
        return gene_symbol_gene_id_and_match_info

    def match_symbols_to_genes(self, release_gene_symbols):
        gene_symbols = (rgs.gene_symbol_id for rgs in release_gene_symbols)
        gene_symbol_gene_id_and_match_info = self._get_gene_id_and_match_info_for_symbol(gene_symbols)

        matches = []
        for release_gene_symbol in release_gene_symbols:
            gene_symbol = release_gene_symbol.gene_symbol_id
            if gene_id_and_match_list := gene_symbol_gene_id_and_match_info.get(gene_symbol):
                for gene_id, match_info in gene_id_and_match_list:
                    matches.append(ReleaseGeneSymbolGene(release_gene_symbol=release_gene_symbol,
                                                         gene_id=gene_id,
                                                         match_info=match_info))

        if matches:
            ReleaseGeneSymbolGene.objects.bulk_create(matches, ignore_conflicts=True, batch_size=2000)

    def match_gene_symbols(self, gene_symbols: Iterable[str]):
        """ gene_symbols must not have been matched  """
        release_gene_symbols = [ReleaseGeneSymbol(release=self.release, gene_symbol_id=gene_symbol_id)
                                for gene_symbol_id in gene_symbols]
        if release_gene_symbols:
            # Need ignore_conflicts=False so we get back PKs
            release_gene_symbols = ReleaseGeneSymbol.objects.bulk_create(release_gene_symbols,
                                                                         batch_size=2000, ignore_conflicts=False)
            self.match_symbols_to_genes(release_gene_symbols)

    def _match_unmatched_gene_symbol_qs(self, gene_symbol_qs):
        """ Match any matched symbols without matched genes """
        unmatched_symbols_qs = gene_symbol_qs.exclude(releasegenesymbol__release=self.release)
        unmatched_symbols = list(unmatched_symbols_qs.values_list("symbol", flat=True).distinct())
        self.match_gene_symbols(unmatched_symbols)

    def match_unmatched_symbols(self, gene_symbol_list):
        gene_symbol_qs = GeneSymbol.objects.filter(pk__in=gene_symbol_list)
        return self._match_unmatched_gene_symbol_qs(gene_symbol_qs)

    def match_unmatched_in_hgnc_and_gene_lists(self):
        """ Make sure symbols have matched genes for each release so they can be used in analyses """
        q_gene_list = Q(genelistgenesymbol__isnull=False)
        q_hgnc = Q(hgnc__isnull=False)
        gene_symbol_qs = GeneSymbol.objects.filter(q_gene_list | q_hgnc)
        return self._match_unmatched_gene_symbol_qs(gene_symbol_qs)


log = logging.getLogger(__name__)

MAX_GENE_SYMBOL_LENGTH = 100


def tokenize_gene_symbols(text) -> TokenizeResult:
    """ Returns TokenizeResult(valid, oversized) where valid is a set of accepted symbols and
        oversized is a list of raw tokens that exceeded MAX_GENE_SYMBOL_LENGTH. Callers should
        surface the oversized list to the user (e.g. via append_gene_list_warning). """
    text = clean_string(text)
    tokens = re.findall(r'[^,;\s]+', text.upper())
    valid, oversized = partition_oversized_names(tokens)
    return TokenizeResult(valid=set(valid), oversized=oversized)


def partition_oversized_names(names: Iterable[str]) -> tuple[list[str], list[str]]:
    """ Splits names into (valid, oversized) by MAX_GENE_SYMBOL_LENGTH. Anything longer cannot
        fit in the genes_genelistgenesymbol unique index (btree row size limit). """
    valid = []
    oversized = []
    for name in names:
        if name is not None and len(name) > MAX_GENE_SYMBOL_LENGTH:
            oversized.append(name)
        else:
            valid.append(name)
    return valid, oversized


def format_oversized_warning(oversized: list[str]) -> str:
    n = len(oversized)
    example = oversized[0]
    preview = example[:40]
    plural = "entry" if n == 1 else "entries"
    return (f"Skipped {n} oversized {plural} (gene symbols must be ≤ {MAX_GENE_SYMBOL_LENGTH} chars). "
            f"First example: '{preview}...' ({len(example)} chars). "
            f"Check the input is a plain gene list.")


def append_gene_list_warning(gene_list, warning: str):
    """ Append a warning to gene_list.error_message and persist it.
        Modifies the in-memory instance so subsequent caller saves preserve the message. """
    if gene_list.error_message:
        gene_list.error_message = f"{gene_list.error_message}\n{warning}"
    else:
        gene_list.error_message = warning
    type(gene_list).objects.filter(pk=gene_list.pk).update(error_message=gene_list.error_message)
