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
