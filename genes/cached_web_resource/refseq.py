import gzip
import logging

import pandas as pd
import requests
from Bio import Entrez

from annotation.models import CachedWebResource
from genes.models import Gene, GeneSymbol, GeneSymbolAlias
from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource
from library.constants import MINUTE_SECS
from library.django_utils import chunked_queryset


def store_refseq_gene_summary_from_web(cached_web_resource: CachedWebResource):
    retrieve_refseq_gene_summaries()

    refseq_genes_qs = Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ)
    num_summaries = refseq_genes_qs.exclude(summary__isnull=True).exclude(summary='').count()
    cached_web_resource.description = f"{num_summaries} RefSeq genes w/summary"
    cached_web_resource.save()


def store_refseq_gene_info_from_web(cached_web_resource: CachedWebResource):
    GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz"
    r = requests.get(GENE_INFO_URL, stream=True, timeout=MINUTE_SECS)
    f = gzip.GzipFile(fileobj=r.raw)
    gene_info_df = pd.read_csv(f, sep='\t')

    known_gene_symbols = GeneSymbol.get_upper_case_lookup()
    known_symbols = set(known_gene_symbols) | set(GeneSymbolAlias.get_upper_case_lookup())

    GeneSymbolAlias.objects.filter(source=GeneSymbolAliasSource.NCBI).delete()
    symbols_synonyms = gene_info_df[["Symbol", "Synonyms"]]
    symbols_synonyms.loc[:, "Symbol"] = symbols_synonyms.loc[:, "Symbol"].str.upper()
    symbols_synonyms.loc[:, "Synonyms"] = symbols_synonyms.loc[:, "Synonyms"].str.upper()

    gene_symbols = []
    gene_symbol_aliases = []
    for _, (symbol, synonyms) in symbols_synonyms.iterrows():
        gene_symbol_id = known_gene_symbols.get(symbol)
        if gene_symbol_id is None:
            gene_symbols.append(GeneSymbol(symbol=symbol))
            gene_symbol_id = symbol

        if synonyms == "-":
            continue

        for s in synonyms.split("|"):
            if s not in known_symbols:
                s = s.upper()
                known_symbols.add(s)
                if s == gene_symbol_id:
                    continue  # Our aliases are case insensitive so no need to store these
                gene_symbol_aliases.append(GeneSymbolAlias(alias=s,
                                                           gene_symbol_id=gene_symbol_id,
                                                           source=GeneSymbolAliasSource.NCBI))

    if gene_symbols:
        GeneSymbol.objects.bulk_create(gene_symbols, ignore_conflicts=True)
    if gene_symbol_aliases:
        GeneSymbolAlias.objects.bulk_create(gene_symbol_aliases, ignore_conflicts=True)

    cached_web_resource.description = f"{len(gene_symbol_aliases)} new gene symbol aliases"
    cached_web_resource.save()


def retrieve_refseq_gene_summaries():
    # 10k limit of return data from NCBI
    BATCH_SIZE = 2000

    refseq_genes_qs = Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ, summary__isnull=True)
    refseq_genes_qs = refseq_genes_qs.exclude(pk__startswith=Gene.FAKE_GENE_ID_PREFIX)

    for genes_qs in chunked_queryset(refseq_genes_qs, BATCH_SIZE):
        gene_ids = genes_qs.values_list("pk", flat=True)
        try:
            gene_annotation = retrieve_entrez_gene_annotation(gene_ids)
        except RuntimeError:
            logging.error("API failure for gene_ids:")
            logging.error(",".join(gene_ids))
            raise

        gene_records = []
        for annot in gene_annotation:
            gene_id = annot.attributes["uid"]
            summary = annot["Summary"]
            gene_records.append(Gene(pk=gene_id, summary=summary))

        if gene_records:
            Gene.objects.bulk_update(gene_records, ["summary"])


def retrieve_entrez_gene_annotation(id_list):
    """ Annotates Entrez Gene IDs (not gene symbol) using Bio.Entrez

        Based on Jinghua (Frank) Feng's code - construct gene reference data """

    request = Entrez.epost("gene", id=",".join(id_list))
    result = Entrez.read(request)

    web_env = result["WebEnv"]
    query_key = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=web_env, query_key=query_key)
    document = Entrez.read(data)

    annotations = document["DocumentSummarySet"]["DocumentSummary"]
    return annotations
