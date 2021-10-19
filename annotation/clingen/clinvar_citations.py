#!/usr/bin/env python3

import logging

import pandas as pd

from annotation.models.models import ClinVarCitation, ClinVarCitationsCollection, Citation
from annotation.models.models_enums import CitationSource
from library.utils import invert_dict

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
    existing_citations = {}
    citation = None
    for citation in Citation.objects.all().order_by("pk"):
        existing_citations[citation.unique_code()] = citation

    if citation:
        max_previously_existing_citation_id = citation.pk  # as qs above is in PK order
    else:
        max_previously_existing_citation_id = 0

    citation_sources = invert_dict(dict(CitationSource.choices))
    new_citations_by_key = {}
    for _, row in df.iterrows():
        #print("row: %s" % row)
        cs = row[CITATION_SOURCE]
        citation_source = citation_sources[cs]

        citation = Citation(citation_source=citation_source,
                            citation_id=row[CITATION_ID])

        key = citation.unique_code()
        if key not in existing_citations:
            new_citations_by_key[key] = citation

    # Insert the new citations
    logging.info("Inserting %d citations", len(new_citations_by_key))
    Citation.objects.bulk_create(new_citations_by_key.values(), batch_size=2000)

    # Update hash
    for citation in Citation.objects.filter(pk__gt=max_previously_existing_citation_id):
        existing_citations[citation.unique_code()] = citation

    # Insert ClinVar citations
    rows = []
    for _, row in df.iterrows():
        cs = row[CITATION_SOURCE]
        citation_source = citation_sources[cs]

        wanted_citation = Citation(citation_source=citation_source,
                                   citation_id=row[CITATION_ID])
        citation = existing_citations[wanted_citation.unique_code()]  # Will die if not there

        cvc = ClinVarCitation(clinvar_citations_collection=clinvar_citations_collection,
                              clinvar_variation_id=row[VARIATION_ID],
                              clinvar_allele_id=row[ALLELE_ID],
                              citation=citation)
        rows.append(cvc)

    num_citations = len(rows)
    logging.info("Read %d records, inserting into DB", num_citations)
    ClinVarCitation.objects.bulk_create(rows, batch_size=2000)

    cached_web_resource.description = f"{num_citations} citations."
    cached_web_resource.save()
