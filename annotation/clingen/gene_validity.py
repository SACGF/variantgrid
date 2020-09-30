from dateutil import parser
from io import StringIO
import logging
import requests

from annotation.models import MonarchDiseaseOntology, DiseaseValidity
from annotation.models.models import GeneDiseaseValidity, GeneDiseaseCurator
from annotation.models.models_enums import ClinGenClassification
from library.utils import invert_dict
import pandas as pd


def store_clingen_gene_validity_curations_from_web(cached_web_resource):
    """ Gets latest copy, maps MONDO -> HPO / OMIM terms """
    CLINGEN_DISEASE_HOST_NAME = "https://search.clinicalgenome.org"
    CLINGEN_DISEASE_CSV_URL = f"{CLINGEN_DISEASE_HOST_NAME}/kb/gene-validity.csv"
    if not MonarchDiseaseOntology.objects.exists():
        msg = "No MonarchDiseaseOntology records, you need to import them - see the annotation page."
        raise ValueError(msg)

    # clingen (grch38) to hg19 genes
    r = requests.get(CLINGEN_DISEASE_CSV_URL)
    df = read_clingen_csv(r.text)

    clingen_disease_validity_lookup = invert_dict(dict(ClinGenClassification.CHOICES))
    genes_set = set()

    clingen_curator = GeneDiseaseCurator.objects.create(name='ClinGen',
                                                        css_class='clingen',
                                                        cached_web_resource=cached_web_resource)

    for _, row in df.iterrows():
        gene_symbol = row['GENE SYMBOL']
        genes_set.add(gene_symbol)
        classification = row['CLASSIFICATION']
        classification = clingen_disease_validity_lookup[classification]
        mondo_id = row['DISEASE ID (MONDO)'].replace("MONDO_", "")
        try:
            mondo = MonarchDiseaseOntology.objects.get(pk=mondo_id)
        except MonarchDiseaseOntology.DoesNotExist:
            mondo = None
            logging.error("Could not load mondo_id '%s' - perhaps update your Monarch Disease Ontology annotations?", mondo_id)

        # Some entries in "CLASSIFICATION DATE" are "false" (WTF?)
        # handle them, but also handle any other non-date that may turn up in the future
        classification_date = None
        classification_date_str = row['CLASSIFICATION DATE']
        if classification_date_str != 'false':
            try:
                classification_date = parser.parse(classification_date_str)
            except Exception as e:
                logging.warning("Couldn't parse classification date '%s': %s", classification_date_str, e)
        disease_validity = DiseaseValidity.objects.create(gene_disease_curator=clingen_curator,
                                                          sop=row['SOP'],
                                                          mondo=mondo,
                                                          classification=classification,
                                                          date=classification_date,
                                                          validity_summary_url=row['ONLINE REPORT'])

        GeneDiseaseValidity.objects.create(disease_validity=disease_validity,
                                           gene_symbol_id=gene_symbol)

    cached_web_resource.description = f"{len(df)} curations for {len(genes_set)} genes"
    cached_web_resource.save()


def read_clingen_csv(clingen_csv_text):
    """ CSV file looks like:

        CLINGEN GENE VALIDITY CURATIONS
        FILE CREATED: 2019-01-04
        WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity
        +++++++++++    ++++++++++++++    +++++++++++++
        GENE SYMBOL    GENE ID (HGNC)    DISEASE LABEL
        +++++++++++    ++++++++++++++    +++++++++++++
        A2ML1    HGNC:23336    Noonan syndrome with multiple lentigines
    """

    HEADER_SEPARATOR = "+++++++++++"
    header_separators_found = 0

    csv_lines = []
    for line in clingen_csv_text.split("\n"):
        if line.startswith(HEADER_SEPARATOR):
            header_separators_found += 1
            continue
        if header_separators_found:
            csv_lines.append(line)

    text = "\n".join(csv_lines)
    f = StringIO(text)
    return pd.read_csv(f, sep=',')
