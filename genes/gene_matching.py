from collections import defaultdict
from typing import Iterable, Dict

from django.db.models import F, Q, Subquery, OuterRef
from django.db.models.functions import Upper
from lazy import lazy
import re

from genes.models import GeneSymbol, GeneSymbolAlias, GeneListGeneSymbol, GeneAnnotationRelease, GeneVersion, \
    ReleaseGeneSymbol, ReleaseGeneSymbolGene, HGNC
from genes.models_enums import HGNCStatus, GeneSymbolAliasSource
from library.cache import timed_cache


class GeneSymbolMatcher:
    @lazy
    def _gene_symbol_lookup(self):
        return GeneSymbol.get_upper_case_lookup()

    @lazy
    def _alias_dict(self):
        return GeneSymbolAlias.get_upper_case_lookup()

    @lazy
    def _release_gene_matchers(self):
        return [GeneMatcher(release) for release in GeneAnnotationRelease.objects.all()]

    def _match_symbols_to_genes_in_releases(self):
        for gm in self._release_gene_matchers:
            gm.match_unmatched_in_hgnc_and_gene_lists()

    def get_gene_symbol_id_and_alias_id(self, original_gene_symbol: str):
        uc_original_gene_symbol = original_gene_symbol.upper()
        gene_symbol_id = self._gene_symbol_lookup.get(uc_original_gene_symbol)
        alias_id = None
        if gene_symbol_id is None:
            if gene_symbol_and_alias := self._alias_dict.get(uc_original_gene_symbol):
                gene_symbol_id, alias_id = gene_symbol_and_alias
        return gene_symbol_id, alias_id

    def get_gene_symbol_id(self, original_gene_symbol: str):
        gene_symbol_id, _alias_id = self.get_gene_symbol_id_and_alias_id(original_gene_symbol)
        return gene_symbol_id

    def create_gene_list_gene_symbols_from_text(self, gene_list, gene_text, save=True):
        gene_names = tokenize_gene_symbols(gene_text)
        return self.create_gene_list_gene_symbols(gene_list, gene_names, save=save)

    def create_gene_list_gene_symbols(self, gene_list, gene_names_list, modification_info=None, save=True):
        """ save=False only returns unsaved objects """
        gene_list_gene_symbols = []
        for original_name in gene_names_list:
            gene_symbol_id, alias_id = self.get_gene_symbol_id_and_alias_id(original_name)
            gene_list_gene_symbols.append(GeneListGeneSymbol(gene_list=gene_list,
                                                             original_name=original_name,
                                                             gene_symbol_id=gene_symbol_id,
                                                             gene_symbol_alias_id=alias_id))

        if save:
            GeneListGeneSymbol.objects.bulk_create(gene_list_gene_symbols, ignore_conflicts=True)
            self._match_symbols_to_genes_in_releases()

        return gene_list_gene_symbols


class HGNCMatcher:
    @staticmethod
    @timed_cache(ttl=3600)
    def instance():
        return HGNCMatcher()

    @lazy
    def _aliases(self) -> Dict:
        alias_qs = GeneSymbolAlias.objects.filter(source=GeneSymbolAliasSource.HGNC)
        alias_qs = alias_qs.annotate(uc_alias=Upper("alias"), uc_symbol=Upper("gene_symbol_id"))
        return dict(alias_qs.values_list("uc_alias", "uc_symbol"))

    @lazy
    def _hgnc_by_uc_gene_symbol(self) -> Dict:
        hgnc_qs = HGNC.objects.filter(status=HGNCStatus.APPROVED)
        return {str(hgnc.gene_symbol_id).upper(): hgnc for hgnc in hgnc_qs}

    def match_hgnc(self, gene_symbol: str):
        gene_symbol = gene_symbol.upper()
        if gs := self._aliases.get(gene_symbol):
            gene_symbol = gs
        return self._hgnc_by_uc_gene_symbol.get(gene_symbol)


class GeneMatcher:
    """ Genes and symbols change over time. This creates DB records connecting them for a particular
        GeneAnnotationRelease.

        To actually retrieve matched genes for a release, use GeneAnnotationRelease.genes_for_symbol(s) """
    def __init__(self, release: GeneAnnotationRelease):
        self.release = release

    @lazy
    def genes(self) -> Dict[str, list]:
        gv_qs = GeneVersion.objects.filter(releasegeneversion__release=self.release).annotate(symbol_upper=Upper("gene_symbol"))
        genes_dict = defaultdict(list)
        for symbol_upper, gene_id in gv_qs.values_list("symbol_upper", "gene__identifier"):
            genes_dict[symbol_upper].append(gene_id)
        return genes_dict

    @lazy
    def aliases_dict(self):
        """ Get symbols from other GeneVersions that match genes from our release """
        qs = GeneVersion.objects.filter(gene__in=self.release.get_genes())
        release_symbol = GeneSymbol.objects.filter(geneversion__gene=OuterRef("gene"),
                                                   geneversion__releasegeneversion__release=self.release)
        diff_version_symbols = qs.annotate(release_symbol=Subquery(release_symbol.values("symbol")[:1])).filter(
            ~Q(gene_symbol=F("release_symbol"))).annotate(symbol_upper=Upper("gene_symbol"))

        values = diff_version_symbols.values_list("symbol_upper", "gene__identifier", "version", "genome_build__name")
        genes_dict = defaultdict(list)
        for symbol_upper, gene_id, version, genome_build in values:
            genes_dict[symbol_upper].append((gene_id, version, genome_build))
        return genes_dict

    def match_gene_symbols(self, gene_symbols: Iterable[str]):
        """ Will fail if ReleaseGeneSymbol already exists """
        release_gene_symbols = []
        gene_symbol_gene_id_and_match_info = defaultdict(list)  # list items = (gene_id, match_info)
        for gene_symbol_id in gene_symbols:
            release_gene_symbols.append(ReleaseGeneSymbol(release=self.release,
                                                          gene_symbol_id=gene_symbol_id))
            gene_name = str(gene_symbol_id).upper()

            if gene_id_list := self.genes.get(gene_name):
                for gene_id in gene_id_list:
                    gene_symbol_gene_id_and_match_info[gene_symbol_id].append((gene_id, None))
            elif alias_list := self.aliases_dict.get(gene_name):
                for gene_id, version, genome_build_name in alias_list:
                    match_info = f"Gene v{version}/{genome_build_name}"
                    gene_symbol_gene_id_and_match_info[gene_symbol_id].append((gene_id, match_info))

        if release_gene_symbols:
            # Need ignore_conflicts=False so we get back PKs
            release_gene_symbols = ReleaseGeneSymbol.objects.bulk_create(release_gene_symbols, batch_size=2000)
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

    def _match_unmatched_gene_symbol_qs(self, gene_symbol_qs):
        """ Match any matched symbols without matched genes """
        unmatched_symbols_qs = gene_symbol_qs.exclude(releasegenesymbol__release=self.release)
        unmatched_symbols = unmatched_symbols_qs.values_list("symbol", flat=True).distinct()
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


def tokenize_gene_symbols(text):
    """ returns set of strings """
    return set(re.findall(r'[^,;\s]+', text.upper()))
