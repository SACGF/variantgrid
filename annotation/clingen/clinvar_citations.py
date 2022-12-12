#!/usr/bin/env python3

import logging
from typing import Dict, Set

import pandas as pd

from annotation.models import Citation
from annotation.models.models import ClinVarCitation, ClinVarCitationsCollection
from annotation.models.models_citations import CitationIdNormalized

ALLELE_ID = '#AlleleID'
VARIATION_ID = 'VariationID'
CITATION_SOURCE = 'citation_source'
CITATION_ID = 'citation_id'
CITATIONS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"


def store_clinvar_citations_from_web(cached_web_resource):
    df = pd.read_csv(CITATIONS_URL, sep='\t', index_col=None)
    for col in [ALLELE_ID, VARIATION_ID, CITATION_SOURCE, CITATION_ID]:
        if col not in df.columns:
            msg = f"Expected column '{col}' in tsv from {CITATIONS_URL}"
            raise ValueError(msg)

    clinvar_citations_collection = ClinVarCitationsCollection.objects.create(cached_web_resource=cached_web_resource)

    rows = list()
    citation_ids: Set[CitationIdNormalized] = set()
    for _, row in df.iterrows():
        citation_source = row[CITATION_SOURCE]

        citation_id = CitationIdNormalized.from_parts(source=citation_source, index=row[CITATION_ID])
        citation_ids.add(citation_id)

        cvc = ClinVarCitation(clinvar_citations_collection=clinvar_citations_collection,
                              clinvar_variation_id=row[VARIATION_ID],
                              clinvar_allele_id=row[ALLELE_ID],
                              citation_id=citation_id.full_id)
        rows.append(cvc)

    logging.info(f"About to ensure existence of {len(citation_ids)} citations")
    Citation.objects.bulk_create(objs=[citation_id.for_bulk_create() for citation_id in citation_ids], batch_size=2000, ignore_conflicts=True)

    num_citations = len(rows)
    logging.info("Read %d records, inserting into DB", num_citations)
    ClinVarCitation.objects.bulk_create(rows, batch_size=2000)

    cached_web_resource.description = f"{num_citations} citations."
    cached_web_resource.save()
