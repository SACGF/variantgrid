"""
`vg inspect gene <symbol|gene id>`: a GeneSymbol (or Gene) - the genes behind the symbol, the latest
GeneVersion per build, aliases in and out, transcript versions per build, canonical transcript
choices, gene lists containing it and classifications naming it.
"""
from typing import Any

from classification.models.classification import Classification
from genes.models.models_gene import Gene, GeneSymbol, GeneVersion, TranscriptVersion
from genes.models.models_gene_coverage import CanonicalTranscript
from genes.models.models_gene_list import GeneListGeneSymbol
from library.vg.inspect import capped, ref
from library.vg.inspect.common import classification_summary
from snpdb.models import GenomeBuild


def load(key: str) -> tuple[GeneSymbol | None, list[Gene]]:
    key = key.strip()
    gene_symbol = GeneSymbol.objects.filter(pk__iexact=key).first()
    if gene_symbol:
        return gene_symbol, list(gene_symbol.get_genes())
    gene = Gene.objects.filter(pk=key).first()
    if gene:
        return None, [gene]
    raise LookupError(f"No GeneSymbol or Gene {key!r}")


def inspect(key: str, depth: int) -> dict[str, Any]:
    gene_symbol, genes = load(key)
    data: dict[str, Any] = {"id": gene_symbol.pk if gene_symbol else genes[0].pk}
    if gene_symbol:
        data["symbol"] = gene_symbol.symbol
        data["url"] = gene_symbol.get_absolute_url()
        meta = gene_symbol.alias_meta
        data["aliases_in"] = [_alias(a) for a in meta.aliases_in][:10]
        data["aliases_out"] = [_alias(a) for a in meta.aliases_out][:10]
    data["genes"] = [_gene(gene) for gene in genes[:10]]
    if depth >= 2:
        transcript_versions = TranscriptVersion.objects.filter(gene_version__gene__in=genes).select_related("transcript", "genome_build")
        data["transcript_versions"] = {
            build.name: capped(transcript_versions.filter(genome_build=build).order_by("transcript_id", "version"),
                               lambda tv: f"{tv.accession}" + (" (MANE)" if getattr(tv, "canonical_score", 0) else ""))
            for build in GenomeBuild.builds_with_annotation().order_by("name")
        }
        if gene_symbol:
            data["canonical_transcripts"] = capped(CanonicalTranscript.objects.filter(gene_symbol=gene_symbol).select_related("collection", "transcript_version__transcript"),
                                                   lambda ct: {"collection": ct.collection.description or ct.collection.filename, "transcript": str(ct.transcript_version) if ct.transcript_version_id else ct.transcript_id})
            data["gene_lists"] = capped(GeneListGeneSymbol.objects.filter(gene_symbol=gene_symbol).select_related("gene_list").order_by("gene_list_id"),
                                        lambda g: ref("gene_list", g.gene_list, g.gene_list.name))
            data["classifications"] = capped(Classification.objects.filter(allele_info__grch38__gene_symbol=gene_symbol).order_by("pk") |
                                             Classification.objects.filter(allele_info__grch37__gene_symbol=gene_symbol).order_by("pk"), classification_summary)
    return data


def _alias(summary) -> str:
    return f"{summary.other_symbol} ({summary.source})" if getattr(summary, "source", None) else str(summary.other_symbol)


def _gene(gene: Gene) -> dict[str, Any]:
    info: dict[str, Any] = {"id": gene.pk, "consortium": gene.get_annotation_consortium_display(), "versions": {}}
    for genome_build in GenomeBuild.builds_with_annotation().order_by("name"):
        version = GeneVersion.objects.filter(gene=gene, genome_build=genome_build).order_by("-version").select_related("gene_symbol").first()
        if version:
            info["versions"][genome_build.name] = {"version": version.version, "symbol": version.gene_symbol_id, "biotype": version.biotype,
                                                   "hgnc": version.hgnc_identifier}
    return info
