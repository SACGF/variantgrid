import gzip

from django.core.management.base import BaseCommand

from genes.cached_web_resource.refseq import store_transcript_sequence_info_fasta
from genes.models_enums import AnnotationConsortium


class Command(BaseCommand):
    """
        Bulk-load TranscriptVersionSequenceInfo from a local FASTA of transcript sequences
        (>accession e.g. NM_000053.3, then the sequence).

        Companion to the "RefSeq Sequence Info" CachedWebResource, which only provides the *latest*
        version of each transcript. Historical classifications often cite older versions - fetching those
        one at a time from NCBI Entrez is slow and rate-limited (see fix_variant_matching --extra), so we
        can persist the fetched stragglers to a supplementary FASTA and reload them here on future deploys.
    """
    def add_arguments(self, parser):
        parser.add_argument('fasta', help='FASTA file (optionally gzipped) of transcript sequences')
        parser.add_argument('--annotation-consortium', default=AnnotationConsortium.REFSEQ.label,
                            choices=[ac.label for ac in AnnotationConsortium])

    def handle(self, *args, **options):
        filename = options['fasta']
        annotation_consortium = {ac.label: ac for ac in AnnotationConsortium}[options['annotation_consortium']]

        handle = gzip.open(filename, "rt") if filename.endswith(".gz") else open(filename, "rt")
        with handle:
            num_records = store_transcript_sequence_info_fasta(handle, filename, annotation_consortium)
        print(f"Loaded {num_records} TranscriptVersionSequenceInfo records from {filename}")
