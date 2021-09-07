"""

Per-gene loss-of-function constraint

SA Path use gnomad_oe_lof EKey which looks like:

0.12(0.06-0.28)
0.03 (0.01 - 0.13)


"""
import gzip
import requests

from genes.models import GnomADGeneConstraint,\
    Transcript, GeneVersion
from genes.models_enums import AnnotationConsortium
import pandas as pd

from library.pandas_utils import df_nan_to_none


def store_gnomad_gene_constraint_from_web(cached_web_resource):
    """ https://gnomad.broadinstitute.org/downloads#gene-constraint """

    GNOMAD_GENE_CONSTRAINT_URL = "https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"

    r = requests.get(GNOMAD_GENE_CONSTRAINT_URL, stream=True)
    f = gzip.GzipFile(fileobj=r.raw)
    df = pd.read_csv(f, sep='\t')
    store_gnomad_gene_constraint_from_df(cached_web_resource, df)


def store_gnomad_gene_constraint_from_df(cached_web_resource, df):
    df = df_nan_to_none(df)

    # gnomAD use Ensembl genes
    transcripts_qs = Transcript.objects.filter(annotation_consortium=AnnotationConsortium.ENSEMBL)
    ensembl_transcript_ids = set(transcripts_qs.values_list("pk", flat=True))

    gene_symbols = set()
    ensembl_gene_ids = set()
    gene_version_qs = GeneVersion.objects.filter(gene__annotation_consortium=AnnotationConsortium.ENSEMBL)

    for gene_id, gene_symbol in gene_version_qs.distinct().values_list("gene_id", "gene_symbol_id"):
        ensembl_gene_ids.add(gene_id)
        gene_symbols.add(gene_symbol)

    gene_constraints = []
    for _, row in df.iterrows():
        gene_id_str = row["gene_id"]
        transcript_id_str = row["transcript"]
        if transcript_id_str in ensembl_transcript_ids:
            transcript_id = transcript_id_str
        else:
            transcript_id = None

        if gene_id_str in ensembl_gene_ids:
            gene_id = gene_id_str
        else:
            gene_id = None

        ggc = GnomADGeneConstraint(gene_symbol_id=row["gene"],
                                   gene_id=gene_id,
                                   transcript_id=transcript_id,
                                   cached_web_resource=cached_web_resource,
                                   oe_lof=row["oe_lof"],
                                   oe_lof_lower=row["oe_lof_lower"],
                                   oe_lof_upper=row["oe_lof_upper"])
        gene_constraints.append(ggc)

    GnomADGeneConstraint.objects.bulk_create(gene_constraints)
    cached_web_resource.description = f"{len(gene_constraints)} genes."
    cached_web_resource.save()
