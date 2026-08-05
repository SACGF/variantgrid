from collections import defaultdict
from typing import Optional

from django.core.management import BaseCommand, CommandError
from django.db.models import Q

from classification.classification_import import reattempt_variant_matching
from classification.models.classification_variant_info_models import ImportedAlleleInfo
from genes.models import TranscriptVersion, TranscriptVersionSequenceInfo
from genes.models_enums import AnnotationConsortium
from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        RefSeq transcripts are sometimes aligned to the genome with gaps. Legacy (pre-cdot) transcript data
        had no gap information, so anything resolved with it could be matched to the wrong variant
        coordinate - cdot data carries the gaps (@see TranscriptVersion.alignment_gap).

        ImportedAlleleInfo caches HGVS resolution, so find the historical ones on gapped transcripts and
        hard rematch them. Selection is deliberately a superset - an ImportedAlleleInfo only records the
        cdot version it was resolved with from #1321 onwards, and re-resolving a record that was already
        right lands on the same variant.
    """

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help='Report the gapped transcripts and counts without rematching')

    @staticmethod
    def _imported_transcript_accession(allele_info: ImportedAlleleInfo) -> Optional[str]:
        if c_hgvs := allele_info.imported_c_hgvs_obj:
            return c_hgvs.transcript
        return allele_info.imported_transcript

    def _gapped_transcript_version(self, genome_build: GenomeBuild, transcript_accession: str) \
            -> Optional[TranscriptVersion]:
        """ The version the record used if we have it, otherwise any version we hold - gaps are a property
            of how the transcript aligns, and versions of a gapped transcript often share the gap """
        try:
            transcript_id, version = TranscriptVersion.get_transcript_id_and_version(transcript_accession)
        except ValueError as ve:
            print(f"Skipping {transcript_accession}: {ve}")
            return None

        tv_qs = TranscriptVersion.objects.filter(transcript_id=transcript_id, genome_build=genome_build)
        transcript_versions = list(tv_qs.order_by("-version"))
        if exact_version := [tv for tv in transcript_versions if tv.version == version]:
            transcript_versions = exact_version

        for transcript_version in transcript_versions:
            try:
                if transcript_version.alignment_gap:
                    return transcript_version
            except Exception as e:
                # Gap detection falls back to transcript length, which needs sequence info from the Entrez API
                print(f"Couldn't determine alignment gap for '{transcript_version.accession}': {e}")
        return None

    def handle(self, *args, **options):
        if not TranscriptVersion.data_is_current_cdot_format():
            raise CommandError(
                "Transcript data is in the old pre-cdot format (TranscriptVersion.data has no "
                "'genome_builds' key), so there is no alignment gap data to compare against. "
                "Run 'python3 manage.py import_cdot_latest' first.")

        # No recorded cdot data version means we can't show it was resolved with gap aware transcript data
        legacy_qs = ImportedAlleleInfo.objects.filter(Q(hgvs_converter_version__isnull=True) |
                                                      Q(hgvs_converter_data_version=""))
        allele_info_ids_by_transcript = defaultdict(list)
        for allele_info in legacy_qs.iterator(chunk_size=1000):
            if genome_build := allele_info.imported_genome_build:
                if transcript_accession := self._imported_transcript_accession(allele_info):
                    allele_info_ids_by_transcript[(genome_build, transcript_accession)].append(allele_info.pk)
        print(f"{len(allele_info_ids_by_transcript)} transcripts to check from legacy ImportedAlleleInfo")

        # Batch retrieve sequence info, otherwise the gap check below falls back to a per-transcript
        # Entrez fetch (2-9s each) for every uncached RefSeq transcript
        refseq_accessions = set()
        for _, transcript_accession in allele_info_ids_by_transcript:
            try:
                if AnnotationConsortium.get_from_transcript_accession(transcript_accession) == AnnotationConsortium.REFSEQ:
                    refseq_accessions.add(transcript_accession)
            except ValueError as ve:
                print(f"Skipping {transcript_accession}: {ve}")

        if refseq_accessions:
            print(f"Batch retrieving sequence info for {len(refseq_accessions)} RefSeq transcripts...")
            TranscriptVersionSequenceInfo.get_refseq_transcript_versions(sorted(refseq_accessions),
                                                                        entrez_batch_size=100,
                                                                        fail_on_error=False)

        allele_info_ids_to_rematch = []
        for (genome_build, transcript_accession), allele_info_ids in allele_info_ids_by_transcript.items():
            if transcript_version := self._gapped_transcript_version(genome_build, transcript_accession):
                print(f"{transcript_accession} ({genome_build}): {transcript_version.cdna_match_diff or 'gap'}"
                      f" - {len(allele_info_ids)} to rematch")
                allele_info_ids_to_rematch.extend(allele_info_ids)

        rematch_qs = ImportedAlleleInfo.objects.filter(pk__in=allele_info_ids_to_rematch)
        print(f"Have {len(allele_info_ids_to_rematch)} ImportedAlleleInfo on gapped transcripts")
        if options["dry_run"]:
            print("Dry run - not rematching")
            return

        reattempt_variant_matching(admin_bot(), rematch_qs, clear_existing=True)
