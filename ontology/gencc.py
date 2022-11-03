from datetime import datetime
from io import StringIO

import requests
from dateutil import parser
from django.conf import settings
from pandas import DataFrame, read_csv

from annotation.models import CachedWebResource
from library.pandas_utils import df_nan_to_none
from library.utils import md5sum_str
from ontology.models import OntologyRelation, GeneDiseaseClassification
from ontology.ontology_builder import OntologyBuilder


def store_gencc_from_web(cached_web_resource: CachedWebResource):
    GENCC_CSV_URL = "https://search.thegencc.org/download/action/submissions-export-csv"
    r = requests.get(GENCC_CSV_URL)
    file_hash = md5sum_str(r.text)
    fileobj = StringIO(r.text)
    num_gene_disease_relationships = load_gencc(fileobj, file_hash, force=True, url=GENCC_CSV_URL)
    cached_web_resource.description = f"{num_gene_disease_relationships} gene/disease relationships"
    cached_web_resource.save()


def _clean_date(date_str) -> datetime.date:
    """ Strip time off date """
    dt = parser.parse(date_str)
    return dt.date()


def load_gencc(file_or_filename, file_hash: str, force: bool, url: str = None) -> int:
    if isinstance(file_or_filename, str):
        url_or_filename = file_or_filename
    else:
        if url is None:
            raise ValueError("If not loading from filename, you must pass in the 'url' parameter")
        url_or_filename = url

    ontology_builder = OntologyBuilder(filename=url_or_filename,
                                       context="gencc_file",
                                       import_source="gencc",
                                       processor_version=7,
                                       force_update=force)
    ontology_builder.ensure_hash_changed(data_hash=file_hash)  # don't re-import if hash hasn't changed
    ontology_builder.cache_everything()

    print("About to load gencc")
    gencc_df: DataFrame = df_nan_to_none(read_csv(file_or_filename, sep=","))

    # only want strong and definitive relationships
    # gencc_df = gencc_df[gencc_df.classification_title.isin(["Definitive", "Strong"])]
    gencc_df = gencc_df.sort_values(by=["gene_curie", "disease_curie"])
    gencc_grouped = gencc_df.groupby(["gene_curie", "disease_curie"])
    gene_disease_classification_lookup = {e.label: e for e in GeneDiseaseClassification}

    num_gene_disease_relationships = 0
    for group_name, df_group in gencc_grouped:
        gene_id = df_group["gene_curie"].iloc[0]
        mondo_id = df_group["disease_curie"].iloc[0]
        gene_symbol = df_group["gene_symbol"].iloc[0]
        sources = []
        classifications = set()

        for row_index, row in df_group.iterrows():
            submitter_title = row["submitter_title"]
            if settings.GENE_RELATION_PANEL_APP_LIVE_UPDATE and submitter_title == "PanelApp Australia":
                continue

            evaluated_date = _clean_date(row["submitted_as_date"])
            submitted_date = _clean_date(row["submitted_run_date"])
            pubmed_ids = []
            if pmids_str := row["submitted_as_pmids"]:
                pubmed_ids = [p.strip() for p in pmids_str.split(",")]

            gencc_classification = row["classification_title"]
            classification_enum = gene_disease_classification_lookup[gencc_classification]
            classifications.add(classification_enum)

            sources.append({
                "submitter": submitter_title,
                "gencc_classification": gencc_classification,
                "mode_of_inheritance":  row["moi_title"],
                "notes": row["submitted_as_notes"],
                "pubmed_ids": pubmed_ids,
                "public_report_url": row["submitted_as_public_report_url"],
                "assertion_criteria_url": row["submitted_as_assertion_criteria_url"],
                "evaluated_date": str(evaluated_date),
                "submitted_date": str(submitted_date),
            })

        if sources:
            num_gene_disease_relationships += 1
            ontology_builder.add_term(
                term_id=gene_id,
                name=gene_symbol,
                definition=None,
                primary_source=False,
                trusted_source=False  # pretty sure GenCC always refers to actual terms, but just in case
            )
            ontology_builder.add_ontology_relation(
                source_term_id=mondo_id,
                dest_term_id=gene_id,
                relation=OntologyRelation.RELATED,
                extra={
                    "sources": sources,
                    "strongest_classification": max(classifications).label,
                }
            )

    ontology_builder.complete(verbose=True)
    return num_gene_disease_relationships
