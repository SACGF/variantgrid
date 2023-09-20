import re
from typing import List

from genes.models import Transcript, TranscriptVersion
from library.utils import first
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

    transcript_version_results: List[TranscriptVersion] = []

    if transcript := Transcript.objects.filter(identifier=transcript_id).first():
        if version:
            transcript_version_for_any_build = False
            for genome_build in search_input.genome_builds:
                if transcript_version := TranscriptVersion.objects.filter(transcript=transcript, version=version, genome_build=genome_build).first():
                    transcript_version_results.append(transcript_version)
                    transcript_version_for_any_build = True

            if not transcript_version_for_any_build:
                yield transcript, f"Unknown transcript version \"{version}\", see transcript page for available versions."
        else:
            yield transcript

    if transcript_version_results:
        if len(transcript_version_results) == 1:
            yield first(transcript_version_results)
        else:
            preview = transcript_version_results[0].preview
            preview.genome_builds = {tv.genome_build for tv in transcript_version_results}
            yield preview
