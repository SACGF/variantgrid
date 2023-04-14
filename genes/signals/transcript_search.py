# from typing import Any, Type
# from django.dispatch import receiver
# from genes.models import Transcript, TranscriptVersion
# from snpdb.search2 import SearchResponseRecordAbstract, SearchResponse, search_signal, SearchInput
# import re
#
#
# TRANSCRIPT_PATTERN = re.compile(r"^(ENST|NM_|NR_|XR_)\d+\.?\d*$")
#
#
# class SearchResponseTranscript(SearchResponseRecordAbstract[Transcript]):
#
#     @classmethod
#     def result_class(cls) -> Type:
#         return Transcript
#
#
# class SearchResponseTranscriptVersion(SearchResponseRecordAbstract[TranscriptVersion]):
#
#     @classmethod
#     def result_class(cls) -> Type:
#         return TranscriptVersion
#
#
# @receiver(search_signal, sender=SearchInput)
# def search_transcript(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
#
#     if search_input.matches_pattern(TRANSCRIPT_PATTERN):
#         response: SearchResponse[SearchResponseTranscript] = SearchResponse(SearchResponseTranscript, search_input)
#         """ return Transcript or TranscriptVersion (build independent) """
#         transcript_id, version = TranscriptVersion.get_transcript_id_and_version(search_input.search_string)
#         obj = None
#         message = None
#         if version:
#             obj = TranscriptVersion.objects.filter(transcript_id=transcript_id, version=version).first()
#             if obj is None:
#                 message = f"Unknown transcript version {version}, see transcript page for available versions."
#
#         if obj is None:  # No version specified or loading version failed
#             obj = Transcript.objects.filter(identifier=transcript_id).first()
#
#         if obj:
#             return [SearchResult(obj, message=message)]
#         return response