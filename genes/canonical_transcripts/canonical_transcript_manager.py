from typing import Set, Tuple

from django.conf import settings

from genes.models import CanonicalTranscriptCollection


class CanonicalTranscriptManager:

    def __init__(self, use_system_default=True):
        """ use_system_default to allow returning None rather than use fallback """
        if use_system_default:
            default_canonical_id = settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID
            if default_canonical_id is None:
                msg = "You need to set settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID"
                raise ValueError(msg)

            try:
                default_canonical_collection = CanonicalTranscriptCollection.objects.get(pk=default_canonical_id)
            except CanonicalTranscriptCollection.DoesNotExist:
                msg = f"Could not load CanonicalTranscriptCollection for settings.GENES_DEFAULT_CANONICAL_TRANSCRIPT_COLLECTION_ID={default_canonical_id}"
                raise ValueError(msg)
        else:
            default_canonical_collection = None
        self.default_canonical_collection = default_canonical_collection
        self.canonical_collection_transcript_lookup = {}

    def get_default_canonical_collection(self):
        return self.default_canonical_collection

    def get_canonical_collection_for_enrichment_kit(self, enrichment_kit):
        canonical_transcript_collection = None
        if enrichment_kit:
            canonical_transcript_collection = enrichment_kit.canonical_transcript_collection

        return canonical_transcript_collection or self.default_canonical_collection

    def get_canonical_transcripts(self, canonical_collection):
        canonical_transcripts = self.canonical_collection_transcript_lookup.get(canonical_collection.pk)
        if canonical_transcripts is None:
            canonical_transcripts = CanonicalTranscriptManager._create_canonical_transcripts(canonical_collection)
            self.canonical_collection_transcript_lookup[canonical_collection.pk] = canonical_transcripts
        return canonical_transcripts

    @staticmethod
    def _create_canonical_transcripts(canonical_collection) -> Tuple[Set, Set]:
        transcript_ids = set()
        original_transcript_ids = set()

        qs = canonical_collection.canonicaltranscript_set.all()
        for transcript_id, original_transcript_id in qs.values_list("transcript_id", "original_transcript_id"):
            transcript_ids.add(transcript_id)
            original_transcript_ids.add(original_transcript_id)
        return transcript_ids, original_transcript_ids
