import logging
from collections import defaultdict

from Bio import SeqIO
from django.core.management import BaseCommand
from django.db.models import QuerySet

from genes.models import TranscriptVersionInfo, TranscriptVersionInfoFastaFileImport, TranscriptVersion, Transcript
from genes.models_enums import AnnotationConsortium
from library.file_utils import file_md5sum, open_handle_gzip




class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--overwrite', action='store_true', help='Delete and replace fasta import with same md5sum')
        parser.add_argument('filename')

    def handle(self, *args, **options):
        filename = options["filename"]
        overwrite = options["overwrite"]

        md5_hash = file_md5sum(filename)
        if existing_import := TranscriptVersionInfoFastaFileImport.objects.filter(md5_hash=md5_hash).first():
            if overwrite:
                print(f"Deleting existing TranscriptVersionInfos for fasta import {md5_hash}")
                existing_import.delete()
            else:
                raise ValueError(f"Fasta import {md5_hash} exists, use --overwrite to delete old data")

        known_transcripts = set(Transcript.objects.all().values_list("identifier", flat=True))
        if not known_transcripts:
            raise ValueError("No transcripts! Insert them first!")

        fasta_import = TranscriptVersionInfoFastaFileImport.objects.create(md5_hash=md5_hash,
                                                                           annotation_consortium=AnnotationConsortium.REFSEQ,
                                                                           filename=filename)
        skipped_transcripts = 0
        records = []
        with open_handle_gzip(filename, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                transcript_id, version = TranscriptVersion.get_transcript_id_and_version(record.id)
                if transcript_id not in known_transcripts:
                    skipped_transcripts += 1
                    continue

                tvi = TranscriptVersionInfo(transcript_id=transcript_id, version=version,
                                            fasta_import=fasta_import,
                                            sequence=str(record.seq), length=len(record.seq))
                records.append(tvi)

        print(f"Skipped {skipped_transcripts} transcripts not in our database")
        if num_records := len(records):
            print(f"Inserting {num_records} TranscriptVersionInfo records")
            TranscriptVersionInfo.objects.bulk_create(records, ignore_conflicts=True, batch_size=2000)

        TranscriptVersionInfo.set_transcript_version_alignment_gap_if_length_different(records)
