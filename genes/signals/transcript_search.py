import re
from genes.models import Transcript, TranscriptVersion
from snpdb.search2 import search_receiver, SearchInputInstance, SearchExample

TRANSCRIPT_PATTERN = re.compile(r"^(ENST|NM_|NR_|XR_)\d+\.?\d*$")


@search_receiver(
    search_type=Transcript,
    pattern=TRANSCRIPT_PATTERN,
    example=SearchExample(
        note="Transcript starting with NM, NR, XR or ENST, version number optional",
        example="NM_001394372"
    )
)
def search_transcript(search_input: SearchInputInstance):
    """ return Transcript or TranscriptVersion (build independent) """
    transcript_id, version = TranscriptVersion.get_transcript_id_and_version(search_input.search_string)

    if transcript := Transcript.objects.filter(identifier=transcript_id).first():
        if version:
            if transcript_version := TranscriptVersion.objects.filter(transcript=transcript, version=version).first():
                yield transcript_version
            else:
                yield transcript, f"Unknown transcript version {version}, see transcript page for available versions."
        else:
            yield transcript
