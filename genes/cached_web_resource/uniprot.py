import ftplib
import gzip
import logging
from collections import defaultdict
from io import BytesIO, TextIOWrapper
from typing import Dict

from Bio import SwissProt
from django.db.models import Q

from annotation.models import CachedWebResource
from genes.models import UniProt

FUNCTION_KEY = "FUNCTION: "
TISSUE_SPEC_KEY = "TISSUE SPECIFICITY: "
PATHWAY_KEY = "PATHWAY: "
PID_KEY = "Pathway_Interaction_DB"
REACTOME_KEY = "Reactome"
JOIN_STRING = " | "  # What's used to merge multiple entries together


def store_uniprot_from_web(cached_web_resource: CachedWebResource):
    logging.info("Downloading 'uniprot_sprot_human.dat.gz' via FTP")
    ftp = ftplib.FTP("ftp.uniprot.org")
    ftp.login("anonymous", "anonymous")
    buffer = BytesIO()
    ftp.retrbinary('RETR /pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_human.dat.gz', buffer.write)
    buffer.seek(0)

    with gzip.GzipFile(fileobj=buffer) as f:
        text_f = TextIOWrapper(f)
        store_uniprot(cached_web_resource, text_f)


def store_uniprot(cached_web_resource: CachedWebResource, file_object):
    logging.info("Extracting data")
    uniprot_data = extract_uniprot_sprot(file_object)

    logging.info("Merging records")
    uniprot_records = []
    for accession, data in uniprot_data.items():
        if data:
            # All values are sets
            str_data = {k: JOIN_STRING.join(v) for k, v in data.items()}
            uniprot_records.append(UniProt(accession=accession, cached_web_resource=cached_web_resource, **str_data))
        else:
            print(f"{accession} had no data we care about!")

    logging.info("Creating DB records")
    UniProt.objects.bulk_create(uniprot_records, batch_size=2000)

    link_hgnc_and_uniprot()

    cached_web_resource.description = f"{len(uniprot_records)} UniProt records"
    cached_web_resource.save()


def extract_uniprot_sprot(f) -> Dict:
    """ Based on Jinghua (Frank) Feng's code - construct gene reference data  """
    reader = SwissProt.parse(f)

    uniprot = defaultdict(lambda: defaultdict(set))
    for record in reader:
        # only use Primary (citable) accession number
        accession = record.accessions[0]
        if accession in uniprot:
            print(f"Adding to existing accession: {accession}")

        for c in record.comments:
            if c.startswith(FUNCTION_KEY):
                uniprot[accession]["function"].add(c.replace(FUNCTION_KEY, '', 1))  # Delete Function_key
            elif c.startswith(TISSUE_SPEC_KEY):
                uniprot[accession]["tissue_specificity"].add(c.replace(TISSUE_SPEC_KEY, '', 1))
            elif c.startswith(PATHWAY_KEY):
                uniprot[accession]["pathway"].add(c.replace(PATHWAY_KEY, '', 1))

        for xfre in record.cross_references:
            if xfre[0] == PID_KEY:
                # ('Pathway_Interaction_DB', 'aurora_a_pathway', 'Aurora A signaling')
                uniprot[accession]["pathway_interaction_db"].add(xfre[2])
            elif xfre[0] == REACTOME_KEY:
                # ('Reactome', 'REACT_111183', 'Meiosis')
                uniprot[accession]["reactome"].add(xfre[2])

    return uniprot


def link_hgnc_and_uniprot():
    """ HGNC links to UniProt objects based on hgnc_ids text field
        the FK is set to NONE on delete (ie during reload of cached_web_resource)

        As HGNC is inserted - it does its own uniprot - so only need to handle case of uniprot being
        inserted last """

    from genes.models import HGNC

    # We only store uni prot IDs that have info we care about, so not all will link
    uniprot_pks = set(UniProt.objects.all().values_list("pk", flat=True))

    hgnc_records = []
    for hgnc in HGNC.objects.filter(uniprot__isnull=True, uniprot_ids__isnull=False).exclude(uniprot_ids=''):
        for up_id in hgnc.uniprot_ids.split(","):
            if up_id in uniprot_pks:
                hgnc.uniprot_id = up_id
                hgnc_records.append(hgnc)

    if hgnc_records:
        logging.info("Linking %d HGNC w/uniprot", len(hgnc_records))
        HGNC.objects.bulk_update(hgnc_records, fields=["uniprot_id"], batch_size=2000)
