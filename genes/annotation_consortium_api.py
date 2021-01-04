from Bio import Entrez
from typing import Tuple
from urllib.error import HTTPError
import requests

from genes.models_enums import AnnotationConsortium
from snpdb.models.models_genome import GenomeBuild


def transcript_exists(genome_build: GenomeBuild, identifier, version) -> Tuple[str, bool]:
    ac_labels = dict(AnnotationConsortium.choices)
    if identifier.startswith("ENST"):
        return ac_labels[AnnotationConsortium.ENSEMBL], ensembl_transcript_exists(genome_build, identifier, version)
    return ac_labels[AnnotationConsortium.REFSEQ], refseq_transcript_exists(identifier, version)


def ensembl_transcript_exists(genome_build: GenomeBuild, identifier, version=None):
    """ If version is passed, check that it's within range of 1 <= version <= latest_version """
    if genome_build.name == 'GRCh37':
        api_base_url = "https://grch37.rest.ensembl.org"
    else:
        api_base_url = "https://rest.ensembl.org"
    url = f"{api_base_url}/lookup/id/{identifier}"
    r = requests.get(url, headers={"Content-Type": "application/json"})
    data = r.json()
    if not r.ok:
        error = data.get("error")
        if error:
            if "not found" in error:
                return False
            raise ValueError(f"Unknown Ensembl API error response: '{error}'")

    if version:
        version = int(version)
        latest_version = int(data["version"])
        return 1 <= version <= latest_version
    return True


def refseq_transcript_exists(identifier, version=None):
    """ No way to tell whether it exists for a particular build or not """
    try:
        accession = identifier
        if version:
            accession += "." + str(version)
        result = Entrez.efetch(db='nuccore', id=accession, retmode='text')
        print(result.read())
        return True
    except HTTPError:
        return False
