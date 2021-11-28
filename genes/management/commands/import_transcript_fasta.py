from collections import Counter

from Bio import SeqIO
from django.core.management import BaseCommand

from genes.models import TranscriptVersionSequenceInfo, TranscriptVersionSequenceInfoFastaFileImport, TranscriptVersion, \
    Transcript
from genes.models_enums import AnnotationConsortium
from library.file_utils import file_md5sum, open_handle_gzip
from library.utils import invert_dict


class Command(BaseCommand):

    def add_arguments(self, parser):
        consortia = [ac[1] for ac in AnnotationConsortium.choices]

        parser.add_argument('--overwrite', action='store_true', help='Delete and replace fasta import with same md5sum')
        parser.add_argument('--annotation-consortium', choices=consortia, required=True)
        parser.add_argument('filename')

    def handle(self, *args, **options):
        filename = options["filename"]
        overwrite = options["overwrite"]
        annotation_consortium_name = options["annotation_consortium"]

        ac_dict = invert_dict(dict(AnnotationConsortium.choices))
        annotation_consortium = ac_dict[annotation_consortium_name]

        md5_hash = file_md5sum(filename)
        if existing_import := TranscriptVersionSequenceInfoFastaFileImport.objects.filter(md5_hash=md5_hash).first():
            if overwrite:
                print(f"Deleting existing TranscriptVersionSequenceInfos for fasta import {md5_hash}")
                existing_import.delete()
            else:
                raise ValueError(f"Fasta import {md5_hash} exists, use --overwrite to delete old data")

        known_transcripts = set(Transcript.objects.all().values_list("identifier", flat=True))
        if not known_transcripts:
            raise ValueError("No transcripts! Insert them first!")

        fasta_import = TranscriptVersionSequenceInfoFastaFileImport.objects.create(md5_hash=md5_hash,
                                                                                   annotation_consortium=annotation_consortium,
                                                                                   filename=filename)
        unknown_transcripts = []
        unknown_transcript_prefixes = Counter()
        records = []
        with open_handle_gzip(filename, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                transcript_id, version = TranscriptVersion.get_transcript_id_and_version(record.id)
                if transcript_id not in known_transcripts:
                    if transcript_id.startswith("X"):
                        continue  # We don't want these
                    prefix = transcript_id.split("_")[0]
                    unknown_transcript_prefixes[prefix] += 1
                    unknown_transcripts.append(Transcript(identifier=transcript_id,
                                                          annotation_consortium=annotation_consortium))

                tvi = TranscriptVersionSequenceInfo(transcript_id=transcript_id, version=version,
                                                    fasta_import=fasta_import,
                                                    sequence=str(record.seq), length=len(record.seq))
                records.append(tvi)

        if unknown_transcripts:
            print(f"Inserting {len(unknown_transcripts)} unknown_transcripts")
            print(unknown_transcript_prefixes)
            Transcript.objects.bulk_create(unknown_transcripts, batch_size=2000)

        if num_records := len(records):
            print(f"Inserting {num_records} TranscriptVersionSequenceInfo records")
            TranscriptVersionSequenceInfo.objects.bulk_create(records, ignore_conflicts=True, batch_size=2000)

        # TranscriptVersionSequenceInfo.set_transcript_version_alignment_gap_if_length_different(records)
