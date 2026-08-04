""" Gene identity for MatchMaker Exchange - one module owning it, so outbound, inbound
    and metrics cannot drift.

    MME's `genomicFeatures[].gene.id` has accepted an Ensembl gene id, an Entrez gene id or
    an HGNC symbol since v1.0; v1.1 puts Ensembl first and makes it mandatory in 2.0. We
    have always emitted a symbol, and compared inbound ids as upper-cased strings - so a
    peer sending ENSG00000133392 for a gene we hold as MYH11 scored zero on the heaviest
    term in the match.

    The fix is not to swap one identifier for another (that just moves which peers we fail),
    but to resolve both sides to a canonical identity and compare on `match_keys()`, while
    publishing Ensembl as `id` alongside the other forms as `_`-prefixed extension fields -
    the worked example in the spec's Non-Standard Fields section.

    Resolution is best-effort by design: annotation completeness varies per deployment and
    MME participation cannot depend on it. Callers treat `None` as an ordinary path.
"""
import logging
import re
import time
from dataclasses import dataclass
from typing import Optional

from django.conf import settings

from classification.enums.classification_enums import SpecialEKeys
from genes.models import GeneSymbol, GeneSymbolAlias, GeneVersion
from genes.models_enums import AnnotationConsortium
from library.constants import HOUR_SECS

# settings.MME_GENE_ID_PREFERENCE - which form goes in the outbound `gene.id`
MME_GENE_ID_ENSEMBL = "ensembl"
MME_GENE_ID_SYMBOL = "symbol"

_ENSEMBL_GENE_PATTERN = re.compile(r"^(ENSG\d+)(?:\.\d+)?$", re.IGNORECASE)
_ENTREZ_GENE_PATTERN = re.compile(r"^(?:GeneID:)?(\d+)$", re.IGNORECASE)


def normalise_gene_key(gene_id: str) -> str:
    """ One spelling per identifier, so `GeneID:4629`/`4629` and `ENSG…`/`ensg….3` compare
        equal. Symbols normalise to upper case, matching ClassificationGroupingSearchTerm. """
    gene_id = gene_id.strip()
    if m := _ENSEMBL_GENE_PATTERN.match(gene_id):
        return m.group(1).upper()
    if m := _ENTREZ_GENE_PATTERN.match(gene_id):
        return m.group(1)
    return gene_id.upper()


@dataclass(frozen=True)
class GeneIdentity:
    """ One gene, in every form MME might name it by. The three fields are what we publish;
        `other_keys` carries the rest of the surface we will match on - superseded symbols
        from older GeneVersions, safe aliases, and any second gene id the symbol maps to. """
    symbol: Optional[str] = None
    ensembl_gene_id: Optional[str] = None
    entrez_gene_id: Optional[str] = None
    other_keys: frozenset[str] = frozenset()

    def match_keys(self) -> set[str]:
        """ Every identifier form this gene is known by, normalised. Both sides of the
            inbound comparison go through this, so the form each peer prefers stops
            mattering. """
        keys = {normalise_gene_key(gene_id) for gene_id
                in (self.symbol, self.ensembl_gene_id, self.entrez_gene_id) if gene_id}
        return keys | set(self.other_keys)

    def canonical_key(self) -> str:
        """ Stable single key for counting - one gene counts once however it was published. """
        return normalise_gene_key(self.ensembl_gene_id or self.symbol or self.entrez_gene_id)

    def as_mme_gene(self) -> dict:
        """ The MME `gene` object. Ensembl as `id` with the alternatives alongside, unless
            the deployment has asked for bare symbols (no Ensembl annotation loaded). """
        if settings.MME_GENE_ID_PREFERENCE != MME_GENE_ID_ENSEMBL:
            return {"id": self.symbol or self.ensembl_gene_id or self.entrez_gene_id}

        gene = {"id": self.ensembl_gene_id or self.symbol or self.entrez_gene_id}
        if self.ensembl_gene_id:
            gene["_ensemblGeneID"] = self.ensembl_gene_id
        if self.entrez_gene_id:
            gene["_entrezGeneID"] = self.entrez_gene_id
        if self.symbol:
            gene["_geneName"] = self.symbol
        return gene


# find_matches() resolves a gene per classification across a scan of up to
# MAX_CLASSIFICATIONS_SCANNED, so resolution is memoised, keyed on the symbol string.
# Gene symbols and aliases only change on a batch import (genes/cached_web_resource/), which
# nothing invalidates caches on - so the whole cache ages out together rather than per
# entry, which is the shape the data actually changes in. Same hour as GeneSymbol._cast.
_SYMBOL_IDENTITY_CACHE: dict[str, Optional[GeneIdentity]] = {}
_symbol_identity_cache_expires: float = 0


def clear_gene_identity_cache():
    """ Also called by tests, which create gene data mid-process. """
    _SYMBOL_IDENTITY_CACHE.clear()


def _identity_cache() -> dict[str, Optional[GeneIdentity]]:
    global _symbol_identity_cache_expires
    if time.time() > _symbol_identity_cache_expires:
        _SYMBOL_IDENTITY_CACHE.clear()
        _symbol_identity_cache_expires = time.time() + HOUR_SECS
    return _SYMBOL_IDENTITY_CACHE


def _identity_for_symbol(gene_symbol: GeneSymbol) -> GeneIdentity:
    """ Add the identifier dimension to a symbol we have already resolved.

        `alias_meta.genes` is the same expansion classification_gene_symbol_filter() uses
        for the datatable - it honours the `different_genes` flag, so aliases that are not
        safe for automatic expansion stay out. Superseded symbols come along because
        GeneVersion rows are per annotation version, which is exactly the inbound case we
        are trying to stop losing: a peer querying by a symbol we have since renamed. """
    cache = _identity_cache()
    cache_key = gene_symbol.symbol.upper()
    if (identity := cache.get(cache_key)) is not None:
        return identity

    genes = {gene for gene in gene_symbol.alias_meta.genes if not gene.is_legacy}
    ensembl_gene_ids = sorted(g.identifier for g in genes
                              if g.annotation_consortium == AnnotationConsortium.ENSEMBL)
    entrez_gene_ids = sorted(g.identifier for g in genes
                             if g.annotation_consortium == AnnotationConsortium.REFSEQ)
    if len(ensembl_gene_ids) > 1:
        logging.debug("MME: gene symbol %s maps to %d Ensembl genes, publishing %s",
                      gene_symbol, len(ensembl_gene_ids), ensembl_gene_ids[0])

    other_keys = set(gene_symbol.alias_meta.alias_symbol_strs)
    other_keys.update(GeneVersion.objects.filter(gene__in=genes, gene_symbol__isnull=False)
                      .values_list("gene_symbol", flat=True).distinct())
    other_keys.update(ensembl_gene_ids)
    other_keys.update(entrez_gene_ids)

    identity = GeneIdentity(
        symbol=gene_symbol.symbol,
        ensembl_gene_id=ensembl_gene_ids[0] if ensembl_gene_ids else None,
        entrez_gene_id=entrez_gene_ids[0] if entrez_gene_ids else None,
        other_keys=frozenset(normalise_gene_key(key) for key in other_keys),
    )
    cache[cache_key] = identity
    return identity


def _identity_for_gene_identifier(identifier: str, annotation_consortium: str) -> GeneIdentity:
    """ Inbound gene id -> the symbol it belongs to, then the full identity. A gene we do
        not hold still self-compares, on the id the peer gave us. """
    gene_version = (GeneVersion.objects
                    .filter(gene__identifier=identifier,
                            gene__annotation_consortium=annotation_consortium,
                            gene_symbol__isnull=False)
                    .select_related("gene_symbol").first())
    if gene_version:
        return _identity_for_symbol(gene_version.gene_symbol)
    if annotation_consortium == AnnotationConsortium.ENSEMBL:
        return GeneIdentity(ensembl_gene_id=identifier)
    return GeneIdentity(entrez_gene_id=identifier)


def _identity_for_symbol_str(symbol_str: str) -> Optional[GeneIdentity]:
    """ GeneSymbol.symbol has a case-insensitive collation, so this matches any spelling.
        A plain query rather than GeneSymbol.cast(), since the identity cache above already
        covers the repeat lookups cast() exists for. """
    cache = _identity_cache()
    cache_key = symbol_str.upper()
    if cache_key in cache:
        return cache[cache_key]

    gene_symbol = GeneSymbol.objects.filter(symbol=symbol_str).first()
    if gene_symbol is None:
        # Symbols outside our GeneVersion history - the alias table is the last route in
        if alias := GeneSymbolAlias.objects.filter(alias=symbol_str).select_related("gene_symbol").first():
            gene_symbol = alias.gene_symbol

    identity = _identity_for_symbol(gene_symbol) if gene_symbol else None
    cache[cache_key] = identity
    return identity


def gene_identity_from_id(gene_id: str) -> Optional[GeneIdentity]:
    """ Resolve any form a peer might send in `genomicFeatures[].gene.id`. """
    gene_id = (gene_id or "").strip()
    if not gene_id:
        return None
    if m := _ENSEMBL_GENE_PATTERN.match(gene_id):
        return _identity_for_gene_identifier(m.group(1).upper(), AnnotationConsortium.ENSEMBL)
    if m := _ENTREZ_GENE_PATTERN.match(gene_id):
        return _identity_for_gene_identifier(m.group(1), AnnotationConsortium.REFSEQ)
    return _identity_for_symbol_str(gene_id)


def resolve_gene(classification) -> Optional[GeneIdentity]:
    """ The gene a classification is about, in every form MME accepts. Starts from the
        curated transcript's resolved symbol where we have one - the same `from_normalized`
        route the grouping search terms use - and falls back to the evidence key. """
    allele_info = getattr(classification, "allele_info", None)
    # the un-normalised per-build FKs rather than resolved_builds - find_matches() scans up
    # to MAX_CLASSIFICATIONS_SCANNED records, and these come along with select_related
    for resolved_variant_info in ((allele_info.grch38, allele_info.grch37) if allele_info else ()):
        if resolved_variant_info and (gene_symbol := resolved_variant_info.gene_symbol):
            return _identity_for_symbol(gene_symbol)
    if symbol_str := classification.get(SpecialEKeys.GENE_SYMBOL):
        return _identity_for_symbol_str(str(symbol_str))
    return None


def gene_match_keys(gene: dict) -> set[str]:
    """ Match surface for one MME `gene` object, resolving `id` and any extension fields a
        peer supplied. Both sides of the inbound Jaccard use this, so a true match produces
        the same key set from either end regardless of which form each published. """
    keys = set()
    for gene_id in (gene.get("id"), gene.get("_ensemblGeneID"),
                    gene.get("_entrezGeneID"), gene.get("_geneName")):
        if not gene_id:
            continue
        if identity := gene_identity_from_id(str(gene_id)):
            keys |= identity.match_keys()
        else:
            keys.add(normalise_gene_key(str(gene_id)))
    return keys


def canonical_gene_key(gene: dict) -> Optional[str]:
    """ One key per gene for counting, so a gene published as ENSG and as a symbol is not
        counted twice. """
    for gene_id in (gene.get("_ensemblGeneID"), gene.get("id"), gene.get("_geneName")):
        if gene_id and (identity := gene_identity_from_id(str(gene_id))):
            return identity.canonical_key()
    if gene_id := gene.get("id"):
        return normalise_gene_key(str(gene_id))
    return None
