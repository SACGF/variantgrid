from urllib.error import HTTPError, URLError

import requests
from Bio import Entrez
from django.core.cache import cache
from requests import RequestException

from genes.models_enums import AnnotationConsortium
from library.constants import WEEK_SECS, HOUR_SECS
from snpdb.models.models_genome import GenomeBuild


def transcript_exists(genome_build: GenomeBuild, identifier, version) -> bool:
    """ Throws exception if we can't tell """

    transcript_exists_key = f"transcript_exists:{identifier}.{version}_{genome_build}"
    transcript_connection_error_key = "error_" + transcript_exists_key
    if cache.get(transcript_connection_error_key):
        raise ConnectionError("Already failed within last hour")

    try:
        # Cache results - true permanently, false for a week and connection error for 1 hour
        exists = cache.get(transcript_exists_key)
        if exists is None:
            if AnnotationConsortium.get_from_transcript_accession(identifier) == AnnotationConsortium.ENSEMBL:
                exists = ensembl_transcript_exists(genome_build, identifier, version)
            else:
                exists = refseq_transcript_exists(identifier, version)
            if exists:
                timeout = None  # Forever
            else:
                timeout = WEEK_SECS
            cache.set(transcript_exists_key, exists, timeout=timeout)
        return exists
    except (RequestException, URLError) as e:
        cache.set(transcript_connection_error_key, True, timeout=HOUR_SECS)
        raise ConnectionError from e


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
        Entrez.efetch(db='nuccore', id=accession, retmode='text')
        return True
    except HTTPError:
        return False
