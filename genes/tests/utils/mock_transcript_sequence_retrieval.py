"""
Serves transcript sequences from recorded FASTAs instead of the RefSeq/Ensembl APIs, so tests never
depend on those being reachable.

To add an accession: fetch it against a real deployment (TranscriptVersionSequenceInfo.get(accession)),
then append '>accession' + its sequence to the FASTA for that consortium.
"""
import os
from typing import Optional

from Bio import SeqIO
from django.conf import settings

from genes.enums import AnnotationConsortium
from genes.transcript_errors import BadTranscript
from genes.transcript_parts import get_transcript_id_and_version
from genes.transcript_sequence_retrieval import FetchedTranscriptSequence, TranscriptSequenceFetcher

TEST_DATA_DIR = os.path.join(settings.BASE_DIR, "genes", "tests", "test_data")
TRANSCRIPT_SEQUENCE_FASTAS = {
    AnnotationConsortium.REFSEQ: "transcript_sequences_refseq.fasta",
    AnnotationConsortium.ENSEMBL: "transcript_sequences_ensembl.fasta",
}


class MockTranscriptSequenceFetcher(TranscriptSequenceFetcher):
    _sequence_by_accession = None

    @classmethod
    def recorded_sequences(cls) -> dict[str, str]:
        if cls._sequence_by_accession is None:
            sequences = {}
            for basename in TRANSCRIPT_SEQUENCE_FASTAS.values():
                with open(os.path.join(TEST_DATA_DIR, basename), encoding="utf-8") as f:
                    for record in SeqIO.parse(f, "fasta"):
                        sequences[record.id] = str(record.seq)
            cls._sequence_by_accession = sequences
        return cls._sequence_by_accession

    @classmethod
    def latest_recorded_version(cls, transcript_id: str) -> Optional[int]:
        versions = []
        for accession in cls.recorded_sequences():
            recorded_id, version = get_transcript_id_and_version(accession)
            if recorded_id == transcript_id and version is not None:
                versions.append(version)
        return max(versions, default=None)

    def fetch(self, transcript_accession: str) -> FetchedTranscriptSequence:
        transcript_id, version = get_transcript_id_and_version(transcript_accession)
        # The APIs resolve an unversioned accession to the latest version
        if version is None:
            version = self.latest_recorded_version(transcript_id)
        accession = f"{transcript_id}.{version}" if version else transcript_accession
        sequence = self.recorded_sequences().get(accession)
        if sequence is None:
            filenames = ", ".join(TRANSCRIPT_SEQUENCE_FASTAS.values())
            raise BadTranscript(f"No recorded sequence for '{transcript_accession}' - add it to one of "
                                f"{filenames} in {TEST_DATA_DIR}")

        return FetchedTranscriptSequence(
            transcript_id=transcript_id,
            annotation_consortium=AnnotationConsortium.get_from_transcript_accession(transcript_accession),
            version=version,
            sequence=sequence,
        )

    def fetch_refseq_batch(self, transcript_accessions, entrez_batch_size: int = 100,
                           fail_on_error=True) -> list[FetchedTranscriptSequence]:
        fetched = []
        for transcript_accession in transcript_accessions:
            try:
                fetched.append(self.fetch(transcript_accession))
            except BadTranscript:
                if fail_on_error:
                    raise
        return fetched
