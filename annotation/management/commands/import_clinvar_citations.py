#!/usr/bin/env python3

from django.core.management.base import BaseCommand
import logging

from annotation.models.models import ClinVarCitation, ClinVarCitationsCollection, \
    Citation
from annotation.models.models_enums import CitationSource
from library.file_utils import file_md5sum
from library.guardian_utils import admin_bot
from library.utils import invert_dict
import pandas as pd
from snpdb.models.models_enums import ImportSource
from upload.models import UploadedFile, UploadedFileTypes, \
    UploadedClinVarCitations

ALLELE_ID = '#AlleleID'
VARIATION_ID = 'VariationID'
CITATION_SOURCE = 'citation_source'
CITATION_ID = 'citation_id'
CITATIONS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/var_citations.txt"


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('var_citations_txt', help=f'File from {CITATIONS_URL}')

    def handle(self, *args, **options):
        filename = options["var_citations_txt"]
        user = admin_bot()

        df = pd.read_csv(filename, sep='\t', index_col=None)
        for col in [ALLELE_ID, VARIATION_ID, CITATION_SOURCE, CITATION_ID]:
            if col not in df.columns:
                msg = f"Expected column '{col}' in tsv from {CITATIONS_URL}"
                raise ValueError(msg)

        logging.info("Deleting existing ClinVarCitations")
        UploadedClinVarCitations.objects.all().delete()

        md5_hash = file_md5sum(filename)
        uploaded_file = UploadedFile.objects.create(path=filename,
                                                    import_source=ImportSource.COMMAND_LINE,
                                                    name='ClinVar citations',
                                                    user=user,
                                                    file_type=UploadedFileTypes.CLINVAR_CITATIONS)

        clinvar_citations_collection = ClinVarCitationsCollection.objects.create()
        UploadedClinVarCitations.objects.create(md5_hash=md5_hash, uploaded_file=uploaded_file, clinvar_citations_collection=clinvar_citations_collection)

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
        Citation.objects.bulk_create(new_citations_by_key.values())

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

        logging.info("Read %d records, inserting into DB", len(rows))
        ClinVarCitation.objects.bulk_create(rows)
