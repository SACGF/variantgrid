import logging

from Bio import Entrez

from annotation.models import CachedWebResource
from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from library.django_utils import chunked_queryset


def store_refseq_gene_summary_from_web(cached_web_resource: CachedWebResource):
    retrieve_refseq_gene_summaries()

    refseq_genes_qs = Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ)
    num_summaries = refseq_genes_qs.exclude(summary__isnull=True).exclude(summary='').count()
    cached_web_resource.description = f"{num_summaries} RefSeq genes w/summary"
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
