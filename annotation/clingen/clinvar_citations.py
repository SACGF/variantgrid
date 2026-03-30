#!/usr/bin/env python3

import io
import logging

import pandas as pd
import requests

from annotation.models import Citation
from annotation.models.models import ClinVarCitation, ClinVarCitationsCollection
from annotation.models.models_citations import CitationIdNormalized

ALLELE_ID = '#AlleleID'
VARIATION_ID = 'VariationID'
CITATION_SOURCE = 'citation_source'
CITATION_ID = 'citation_id'
CITATIONS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"
DOWNLOAD_TIMEOUT_SECONDS = 120


def store_clinvar_citations_from_web(cached_web_resource):
    response = requests.get(CITATIONS_URL, timeout=DOWNLOAD_TIMEOUT_SECONDS)
    response.raise_for_status()

    df = pd.read_csv(io.BytesIO(response.content), sep='\t', index_col=None,
                     dtype={ALLELE_ID: int,
                            VARIATION_ID: int,
                            "rs": str,
                            "nsv": str,
                            CITATION_SOURCE: str,
                            CITATION_ID: str})
    for col in [ALLELE_ID, VARIATION_ID, CITATION_SOURCE, CITATION_ID]:
        if col not in df.columns:
            msg = f"Expected column '{col}' in tsv from {CITATIONS_URL}"
            raise ValueError(msg)

    clinvar_citations_collection = ClinVarCitationsCollection.objects.create(cached_web_resource=cached_web_resource)

    rows = []
    citation_ids: set[CitationIdNormalized] = set()
    skipped = 0
    for _, row in df.iterrows():
        citation_source = row[CITATION_SOURCE]
        citation_id_val = row[CITATION_ID]

        if not citation_source or not citation_id_val or pd.isna(citation_source) or pd.isna(citation_id_val):
            skipped += 1
            continue

        try:
            citation_id = CitationIdNormalized.from_parts(source=citation_source, index=citation_id_val)
        except ValueError:
            logging.warning("Skipping invalid citation row: source=%r id=%r", citation_source, citation_id_val)
            skipped += 1
            continue

        citation_ids.add(citation_id)

        cvc = ClinVarCitation(clinvar_citations_collection=clinvar_citations_collection,
                              clinvar_variation_id=row[VARIATION_ID],
                              clinvar_allele_id=row[ALLELE_ID],
                              citation_id=citation_id.full_id)
        rows.append(cvc)

    if skipped:
        logging.warning("Skipped %d invalid rows from %s", skipped, CITATIONS_URL)

    logging.info("About to ensure existence of %d citations", len(citation_ids))
    Citation.objects.bulk_create(objs=[citation_id.for_bulk_create() for citation_id in citation_ids], batch_size=2000, ignore_conflicts=True)

    num_citations = len(rows)
    logging.info("Read %d records, inserting into DB", num_citations)
    ClinVarCitation.objects.bulk_create(rows, batch_size=2000)

    cached_web_resource.description = f"{num_citations} citations."
    cached_web_resource.save()
