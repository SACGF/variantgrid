from typing import Any, Type
from django.dispatch import receiver
from genes.models import Transcript, TranscriptVersion
from snpdb.search2 import SearchResponse, search_signal, SearchInput
import re


TRANSCRIPT_PATTERN = re.compile(r"^(ENST|NM_|NR_|XR_)\d+\.?\d*$")



@receiver(search_signal, sender=SearchInput)
def search_transcript(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:

    if search_input.matches_pattern(TRANSCRIPT_PATTERN):
        response = SearchResponse(Transcript)
        """ return Transcript or TranscriptVersion (build independent) """
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(search_input.search_string)

        if transcript := Transcript.objects.filter(identifier=transcript_id).first():
            if version:
                if transcript_version := TranscriptVersion.objects.filter(transcript=transcript, version=version).first():
                    response.add(transcript_version)
                else:
                    response.add(transcript, messages=[f"Unknown transcript version {version}, see transcript page for available versions."])
            else:
                response.add(transcript)
        return response