import logging
from datetime import datetime
from io import StringIO

import requests
from dateutil import parser
from django.conf import settings
from pandas import DataFrame, isna, read_csv

from annotation.models import CachedWebResource
from library.constants import MINUTE_SECS
from library.pandas_utils import df_nan_to_none
from library.utils import md5sum_str
from ontology.models import OntologyRelation, GeneDiseaseClassification
from ontology.ontology_builder import OntologyBuilder


def store_gencc_from_web(cached_web_resource: CachedWebResource):
    GENCC_CSV_URL = "https://search.thegencc.org/download/action/submissions-export-csv"
    r = requests.get(GENCC_CSV_URL, timeout=MINUTE_SECS)
    file_hash = md5sum_str(r.text)
    fileobj = StringIO(r.text)
    num_gene_disease_relationships = load_gencc(fileobj, file_hash, force=True, url=GENCC_CSV_URL)
    cached_web_resource.description = f"{num_gene_disease_relationships} gene/disease relationships"
    cached_web_resource.save()


def _clean_date(date_str) -> datetime.date:
    """ Strip time off date """
    dt = parser.parse(date_str)
    return dt.date()


def _none_if_na(value):
    """ pandas iterrows() coerces df_nan_to_none's Python None back to float NaN — re-normalise here. """
    if value is None or isna(value):
        return None
    return value


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

    logging.info("About to load gencc")
    gencc_df: DataFrame = df_nan_to_none(read_csv(file_or_filename, sep=",", dtype={"submitted_as_pmids": str}))

    # only want strong and definitive relationships
    gencc_df = gencc_df.sort_values(by=["gene_curie", "disease_curie"])
    gene_disease_classification_lookup = {e.label: e for e in GeneDiseaseClassification}

    num_gene_disease_relationships = 0
    for _, df_group in gencc_df.groupby(["gene_curie", "disease_curie"]):
        gene_id = df_group["gene_curie"].iloc[0]
        mondo_id = df_group["disease_curie"].iloc[0]
        gene_symbol = df_group["gene_symbol"].iloc[0]
        sources = []
        classifications = set()

        for _, row in df_group.iterrows():
            submitter_title = row["submitter_title"]
            if settings.GENE_RELATION_PANEL_APP_LIVE_UPDATE and submitter_title == "PanelApp Australia":
                continue

            evaluated_date = _clean_date(row["submitted_as_date"])
            submitted_date = _clean_date(row["submitted_run_date"])
            pubmed_ids = []
            pmids_value = row["submitted_as_pmids"]
            if pmids_value is not None and not isna(pmids_value):
                pmids_str = str(pmids_value)
                if pmids_str.endswith(".0"):
                    pmids_str = pmids_str[:-2]
                pubmed_ids = [p.strip() for p in pmids_str.split(",") if p.strip()]

            gencc_classification = row["classification_title"]
            classification_enum = gene_disease_classification_lookup[gencc_classification]
            classifications.add(classification_enum)

            sources.append({
                "submitter": _none_if_na(submitter_title),
                "gencc_classification": gencc_classification,
                "mode_of_inheritance": _none_if_na(row["moi_title"]),
                "notes": _none_if_na(row["submitted_as_notes"]),
                "pubmed_ids": pubmed_ids,
                "public_report_url": _none_if_na(row["submitted_as_public_report_url"]),
                "assertion_criteria_url": _none_if_na(row["submitted_as_assertion_criteria_url"]),
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
