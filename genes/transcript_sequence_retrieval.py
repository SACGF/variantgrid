"""
Retrieval of transcript sequences from the RefSeq (NCBI Entrez) and Ensembl REST APIs.

Fetching only - the caller stores what comes back, so this stays free of the ORM and can be swapped for
an implementation that serves recorded sequences (@see genes.tests.utils.mock_transcript_sequence_retrieval)
"""
import logging
from dataclasses import dataclass
from io import StringIO
from typing import Iterable, Optional
from urllib.error import HTTPError

import requests
from Bio import Entrez, SeqIO

from genes.models_enums import AnnotationConsortium
from genes.transcript_errors import BadTranscript, NoTranscript
from genes.transcript_parts import get_transcript_id_and_version
from library.constants import MINUTE_SECS
from library.utils import get_single_element, iter_fixed_chunks
from snpdb.models import GenomeBuild

ENSEMBL_REST_BASE_URLS = {
    "GRCh37": "https://grch37.rest.ensembl.org",
    "GRCh38": "https://rest.ensembl.org",
}


@dataclass(frozen=True)
class FetchedTranscriptSequence:
    """ A sequence retrieved from an external API, before we store it """
    transcript_id: str
    annotation_consortium: str
    version: int
    sequence: str
    api_response: Optional[str] = None
    # Ensembl only serves the latest version, which may not be the one asked for - we store both
    also_store_version: Optional[int] = None

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def versions(self) -> list[int]:
        versions = [self.version]
        if self.also_store_version is not None and self.also_store_version != self.version:
            versions.append(self.also_store_version)
        return versions


class TranscriptSequenceFetcher:
    """ Retrieves transcript sequences we don't hold locally """

    override_class = None  # Tests sub in a recorded-data implementation - @see variantgrid.test_runner

    @classmethod
    def instance(cls) -> 'TranscriptSequenceFetcher':
        return (cls.override_class or cls)()

    def fetch(self, transcript_accession: str) -> FetchedTranscriptSequence:
        annotation_consortium = AnnotationConsortium.get_from_transcript_accession(transcript_accession)
        if annotation_consortium == AnnotationConsortium.REFSEQ:
            return self.fetch_refseq(transcript_accession)
        return self.fetch_ensembl(transcript_accession)

    def fetch_refseq(self, transcript_accession: str) -> FetchedTranscriptSequence:
        try:
            data = Entrez.efetch(db='nuccore', id=transcript_accession, rettype='gb', retmode='text')
        except HTTPError as e:
            if e.code == 400:
                raise BadTranscript(f"Bad Transcript: Entrez API reports \"{transcript_accession}\" not found") from e
            raise e
        api_response = data.read()
        with StringIO(api_response) as f:
            if api_response.startswith("Error:"):
                error_message = api_response[6:]
                if len(error_message) > 2 and error_message[2] == " ":
                    # error messages seem to have been iterated so that every 2nd characer is a space, fix this
                    error_message = "".join(char for i, char in enumerate(error_message) if i % 2 == 1)
                raise ValueError("ClinGen API response: " + error_message)

            record = get_single_element(list(SeqIO.parse(f, "genbank")))
            return _fetched_from_genbank_record(record, api_response)

    def fetch_ensembl(self, transcript_accession: str,
                      genome_build: GenomeBuild = None) -> FetchedTranscriptSequence:
        if genome_build is None:
            genome_build = GenomeBuild.grch38()
        base_url = ENSEMBL_REST_BASE_URLS[genome_build.name]
        transcript_id, requested_version = get_transcript_id_and_version(transcript_accession)
        url = f"{base_url}/sequence/id/{transcript_id}?type=cdna"
        r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=MINUTE_SECS)
        data = r.json()

        if r.ok:
            return FetchedTranscriptSequence(transcript_id=data["id"],
                                             annotation_consortium=AnnotationConsortium.ENSEMBL,
                                             version=requested_version,
                                             sequence=data["seq"],
                                             api_response=r.text,
                                             also_store_version=data["version"])

        error = data.get("error")
        if error:
            if "not found" in error:
                if genome_build != GenomeBuild.grch37():
                    # Try 37 this time
                    return self.fetch_ensembl(transcript_accession, GenomeBuild.grch37())
                raise BadTranscript(f"Ensembl API reports '{transcript_id}' not found")
        raise NoTranscript(f"Unable to understand Ensembl API response: {data}")

    def fetch_refseq_batch(self, transcript_accessions: Iterable[str], entrez_batch_size: int = 100,
                           fail_on_error=True) -> list[FetchedTranscriptSequence]:
        fetched = []
        for id_list in iter_fixed_chunks(transcript_accessions, entrez_batch_size):
            id_param = ",".join(id_list)
            try:
                search_results = Entrez.read(Entrez.epost("nuccore", id=id_param))
                fetch_handle = Entrez.efetch(
                    db="nuccore",
                    rettype="gb",
                    retmode="text",
                    webenv=search_results["WebEnv"],
                    query_key=search_results["QueryKey"],
                    idtype="acc",
                )
                for record in SeqIO.parse(fetch_handle, "genbank"):
                    # Store raw data so that we can retrieve more stuff from it later
                    s = StringIO()
                    SeqIO.write(record, s, "genbank")
                    s.seek(0)
                    fetched.append(_fetched_from_genbank_record(record, s.read()))
            except Exception as e:
                # Entrez surfaces failures as RuntimeError, urllib HTTPError (e.g. 429 rate limit), etc.
                # When the caller passed fail_on_error=False they want a best-effort batch, so swallow them all.
                logging.warning("Entrez failed w/params: %s (%s)", id_param, e)
                if fail_on_error:
                    raise e
        return fetched


def _fetched_from_genbank_record(record, api_response: str) -> FetchedTranscriptSequence:
    transcript_id, version = get_transcript_id_and_version(record.id)
    return FetchedTranscriptSequence(transcript_id=transcript_id,
                                     annotation_consortium=AnnotationConsortium.REFSEQ,
                                     version=version,
                                     sequence=str(record.seq),
                                     api_response=api_response)
