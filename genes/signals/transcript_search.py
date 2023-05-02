import re
from genes.models import Transcript, TranscriptVersion
from snpdb.search import search_receiver, SearchInputInstance, SearchExample

TRANSCRIPT_PATTERN = re.compile(r"^(ENST|NM_|NR_|XR_)\d+\.?\d*$", re.IGNORECASE)


@search_receiver(
    search_type=Transcript,
    pattern=TRANSCRIPT_PATTERN,
    example=SearchExample(
        note="Transcript starting with NM, NR, XR or ENST, version number optional",
        examples=["NM_001394372"]
    )
)
def search_transcript(search_input: SearchInputInstance):
    """ return Transcript or TranscriptVersion (build independent) """
    upper_string = search_input.search_string.upper()
    transcript_id, version = TranscriptVersion.get_transcript_id_and_version(upper_string)

    if transcript := Transcript.objects.filter(identifier=transcript_id).first():
        if version:
            if transcript_version := TranscriptVersion.objects.filter(transcript=transcript, version=version).first():
                yield transcript_version
            else:
                yield transcript, f"Unknown transcript version {version}, see transcript page for available versions."
        else:
            yield transcript
