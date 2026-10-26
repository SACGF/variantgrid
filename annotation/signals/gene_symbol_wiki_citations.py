"""
gene_symbol_wiki_post_save_handler keeps a gene symbol's GeneSymbolCitations in step with the citations in its
wiki markdown (PubMed, PubMed Central and NCBI Bookshelf references), adding new ones and dropping those no longer
mentioned. Connected in annotation/apps.py.
"""
import html
import logging

from django.utils.html import strip_tags
from markdown import markdown

from annotation.models.models_citations import CitationFetchRequest
from annotation.regexes import db_citation_regexes
from genes.models import GeneSymbolWiki


def gene_symbol_wiki_post_save_handler(sender, instance: GeneSymbolWiki, **kwargs):  # pylint: disable=unused-argument
    gene_symbol = instance.gene_symbol
    existing_citations = {gc.citation_id: gc for gc in gene_symbol.genesymbolcitation_set.all()}

    gene_citation_ids_to_keep = []
    # Text pasted from PubMed arrives formatted ("PMID: **21892162**"), which hides the id from the regex
    plain_text = html.unescape(strip_tags(markdown(instance.markdown)))
    citation_refs = db_citation_regexes.search(plain_text)
    for citation in CitationFetchRequest.get_unfetched_citations(citation_refs):
        gene_citation = existing_citations.get(citation.pk)
        if not gene_citation:
            gene_citation = gene_symbol.genesymbolcitation_set.create(citation=citation)
        gene_citation_ids_to_keep.append(gene_citation.pk)

    num_deleted = gene_symbol.genesymbolcitation_set.exclude(pk__in=gene_citation_ids_to_keep).delete()[0]
    logging.info("Kept %s, deleted %d", gene_citation_ids_to_keep, num_deleted)
