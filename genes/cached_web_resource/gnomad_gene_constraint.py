"""

Per-gene loss-of-function constraint

SA Path use gnomad_oe_lof EKey which looks like:

0.12(0.06-0.28)
0.03 (0.01 - 0.13)


"""
import pandas as pd

from genes.models import GnomADGeneConstraint, GeneSymbol, TranscriptVersion
from library.pandas_utils import df_nan_to_none
from snpdb.models import GenomeBuild


def store_gnomad_gene_constraint_from_web(cached_web_resource):
    """ https://gnomad.broadinstitute.org/downloads#v4-constraint """

    GNOMAD4_GENE_CONSTRAINT_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/v4.0/constraint/gnomad.v4.0.constraint_metrics.tsv"
    df = pd.read_csv(GNOMAD4_GENE_CONSTRAINT_URL, sep='\t')
    store_gnomad_gene_constraint_from_df(cached_web_resource, df)


def store_gnomad_gene_constraint_from_df(cached_web_resource, df):
    df = df_nan_to_none(df)

    # gnomAD4 uses both RefSeq and Ensembl, but Ensembl transcripts don't have versions while RefSeq does
    transcript_versions_by_id = TranscriptVersion.transcript_versions_by_id(GenomeBuild.grch38())

    skipped_records = 0
    gene_symbol_uc_lookup = GeneSymbol.get_upper_case_lookup()

    new_gene_symbols = set()
    gene_constraints = []
    for _, row in df.iterrows():
        original_transcript = row["transcript"]
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(original_transcript)
        transcript = None
        transcript_version_id = None
        if transcript_versions := transcript_versions_by_id.get(transcript_id):
            transcript_version_id = transcript_versions.get(version)
            transcript = transcript_id  # Matched

        gene_symbol = row["gene"]
        if gene_symbol is None or gene_symbol == 'NA':
            skipped_records += 1
            continue

        gene_symbol_id = gene_symbol_uc_lookup.get(gene_symbol.upper())
        if gene_symbol_id is None:
            new_gene_symbols.add(GeneSymbol(symbol=gene_symbol))

        # Convert from gnomAD to django fields is
        skip_fields = {"gene", "transcript"}
        django_fields = {}
        gnomad_field: str
        for gnomad_field, value in row.items():
            if gnomad_field in skip_fields:
                continue
            f = gnomad_field.replace(".", "_").lower()
            django_fields[f] = value

        ggc = GnomADGeneConstraint(gene_symbol_id=gene_symbol,
                                   transcript_id=transcript,
                                   transcript_version_id=transcript_version_id,
                                   cached_web_resource=cached_web_resource,
                                   **django_fields)
        gene_constraints.append(ggc)

    if new_gene_symbols:
        print(f"Inserting {len(new_gene_symbols)} new gene symbols")
        GeneSymbol.objects.bulk_create(new_gene_symbols, batch_size=2000, ignore_conflicts=True)

    if skipped_records:
        print(f"Skipped {skipped_records} records (due to missing gene symbol)")
    GnomADGeneConstraint.objects.bulk_create(gene_constraints)
    cached_web_resource.description = f"{len(gene_constraints)} genes."
    cached_web_resource.save()
