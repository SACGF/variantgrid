from collections import defaultdict
from typing import Optional

from django.core.management import BaseCommand, CommandError

from classification.classification_import import reattempt_variant_matching
from classification.models.classification_variant_info_models import ImportedAlleleInfo
from genes.hgvs import HGVSDisplay
from genes.models import TranscriptVersion, TranscriptVersionSequenceInfo
from genes.models_enums import AnnotationConsortium
from library.guardian_utils import admin_bot
from snpdb.models import GenomeBuild


class Command(BaseCommand):
    """
        RefSeq transcripts are sometimes aligned to the genome with gaps. Legacy (pre-cdot) transcript data
        had no gap information, so anything resolved with it could be matched to the wrong variant
        coordinate - cdot data carries the gaps (@see TranscriptVersion.alignment_gap).

        ImportedAlleleInfo caches HGVS resolution, so find the ones on gapped transcripts whose matched
        variant disagrees with what current transcript data resolves to, and hard rematch those.

        The converter version stamped on an ImportedAlleleInfo isn't evidence of how the variant was
        matched - update_variant_coordinate/recalc_c_hgvs re-stamp it while leaving the matched variant
        as it was, so the stored coordinate is compared against a fresh resolution instead.
    """

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true',
                            help='Report the gapped transcripts and differences without rematching')

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

    @staticmethod
    def _matching_out_of_date(allele_info: ImportedAlleleInfo) -> Optional[str]:
        """ How the matched variant differs from what current transcript data resolves to, None if they agree """
        genome_build = allele_info.imported_genome_build
        new_coordinate = allele_info.calculate_variant_coordinate().variant_coordinate
        if new_coordinate is None:
            return None  # Can't resolve it now, so a rematch could only lose the match we have

        if matched_variant := allele_info.matched_variant:
            old_coordinate = matched_variant.coordinate
            if old_coordinate.as_internal_symbolic(genome_build) != new_coordinate.as_internal_symbolic(genome_build):
                return f"{old_coordinate} -> {new_coordinate}"
            return None
        return f"unmatched -> {new_coordinate}"

    def handle(self, *args, **options):
        if not TranscriptVersion.data_is_current_cdot_format():
            raise CommandError(
                "Transcript data is in the old pre-cdot format (TranscriptVersion.data has no "
                "'genome_builds' key), so there is no alignment gap data to compare against. "
                "Run 'python3 manage.py import_cdot_latest' first.")

        genome_builds_by_name = {gb.name: gb for gb in GenomeBuild.builds_with_annotation()}
        allele_info_ids_by_transcript = defaultdict(list)
        allele_info_values = ImportedAlleleInfo.objects.values_list(
            "pk", "imported_c_hgvs", "imported_transcript", "imported_genome_build_patch_version__genome_build_id")
        for pk, imported_c_hgvs, imported_transcript, genome_build_name in allele_info_values.iterator(chunk_size=10000):
            if genome_build := genome_builds_by_name.get(genome_build_name):
                transcript_accession = HGVSDisplay(imported_c_hgvs).transcript if imported_c_hgvs else imported_transcript
                if transcript_accession:
                    allele_info_ids_by_transcript[(genome_build, transcript_accession)].append(pk)
        print(f"{len(allele_info_ids_by_transcript)} transcripts used by ImportedAlleleInfo to check")

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

        gapped_allele_info_ids = []
        for (genome_build, transcript_accession), allele_info_ids in allele_info_ids_by_transcript.items():
            if transcript_version := self._gapped_transcript_version(genome_build, transcript_accession):
                print(f"{transcript_accession} ({genome_build}): {transcript_version.cdna_match_diff or 'gap'}"
                      f" - {len(allele_info_ids)} ImportedAlleleInfo")
                gapped_allele_info_ids.extend(allele_info_ids)
        print(f"{len(gapped_allele_info_ids)} ImportedAlleleInfo on gapped transcripts - checking their matching")

        allele_info_ids_to_rematch = []
        for allele_info in ImportedAlleleInfo.objects.filter(pk__in=gapped_allele_info_ids).iterator(chunk_size=100):
            if difference := self._matching_out_of_date(allele_info):
                print(f"REMATCHING {allele_info.pk} ({allele_info.imported_hgvs}): "
                      f"{difference}")
                allele_info_ids_to_rematch.append(allele_info.pk)

        rematch_qs = ImportedAlleleInfo.objects.filter(pk__in=allele_info_ids_to_rematch)
        print(f"Have {len(allele_info_ids_to_rematch)} ImportedAlleleInfo to rematch")
        if options["dry_run"]:
            print("Dry run - not rematching")
            return

        reattempt_variant_matching(admin_bot(), rematch_qs, clear_existing=True)
