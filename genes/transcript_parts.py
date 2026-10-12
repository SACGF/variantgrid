from dataclasses import dataclass
from typing import Optional

from library.utils import FormerTuple


@dataclass
class TranscriptParts(FormerTuple):
    identifier: str
    version: Optional[int]

    @property
    def as_tuple(self) -> tuple:
        return self.identifier, self.version

    def __repr__(self):
        if self.version:
            return f"{self.identifier}.{self.version}"
        return self.identifier


def get_transcript_id_and_version(transcript_accession: str) -> TranscriptParts:
    """ Lenient - anything that isn't "<identifier>.<int>" is kept whole as the identifier """
    parts = transcript_accession.split(".")
    if len(parts) == 2 and parts[1].isdigit():
        identifier = str(parts[0])
        version = int(parts[1])
    else:
        identifier, version = transcript_accession, None
    return TranscriptParts(identifier, version)


# RefSeq gives the mitochondrial coding genes no RNA accession, so cdot names them 'fake-rna-<gene>' (eg
# 'fake-rna-ND4') with no version. VEP's RefSeq cache names the same GFF records after the gene ('ND4.1')
CDOT_FAKE_TRANSCRIPT_PREFIX = "fake-rna-"
CDOT_FAKE_TRANSCRIPT_VERSION = 1


def get_cdot_transcript_id_and_version(transcript_accession: str) -> TranscriptParts:
    """ As get_transcript_id_and_version, but a cdot fake transcript gets CDOT_FAKE_TRANSCRIPT_VERSION so it
        can be stored as a TranscriptVersion (#2139) """
    transcript_parts = get_transcript_id_and_version(transcript_accession)
    if transcript_parts.version is None and transcript_accession.startswith(CDOT_FAKE_TRANSCRIPT_PREFIX):
        transcript_parts = TranscriptParts(transcript_accession, CDOT_FAKE_TRANSCRIPT_VERSION)
    return transcript_parts


def get_cdot_fake_transcript_id(gene_name: str) -> str:
    """ The cdot fake transcript for a gene with no RefSeq RNA accession, ie 'ND4' -> 'fake-rna-ND4' """
    return f"{CDOT_FAKE_TRANSCRIPT_PREFIX}{gene_name}"
