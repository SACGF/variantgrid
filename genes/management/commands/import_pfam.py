#!/usr/bin/env python3
# Data can be downloaded here:
# * ftp://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/idmapping/by_organism/
# * ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes
import re

from django.core.management.base import BaseCommand
from genes.models import TranscriptVersion, Pfam, Transcript, PfamSequence, \
    PfamSequenceIdentifier, PfamDomains
from genes.models_enums import AnnotationConsortium
from library.log_utils import console_logger
import pandas as pd

class Command(BaseCommand):
    BULK_INSERT_SIZE = 2000
    PFAM_SEQUENCE_PATTERN = re.compile(r"(.+?)-\d+$")

    def add_arguments(self, parser):
        parser.add_argument('--idmapping', required=True, help="eg HUMAN_9606_idmapping.dat.gz")
        parser.add_argument('--pfam-tsv', required=True, help="eg 9606.tsv.gz")

    def handle(self, *args, **options):
        idmapping = options["idmapping"]
        pfam_tsv = options["pfam_tsv"]
        logger = console_logger()

        if not Pfam.objects.exists():
            raise ValueError("No Pfam in the database! You need to retrieve that from annotation page first!")

        # Delete existing
        PfamSequenceIdentifier.objects.all().delete()
        PfamDomains.objects.all().delete()

        mapping_df = pd.read_csv(idmapping, header=None, names=["seq_id", "id_type", "identifier"], sep='\t')
        ensembl_transcript_mask = mapping_df["id_type"] == "Ensembl_TRS"
        refseq_transcript_mask = mapping_df["id_type"] == "RefSeq_NT"

        # Insert sequences
        unique_sequences = mapping_df[ensembl_transcript_mask | refseq_transcript_mask]["seq_id"].unique()
        self.insert_sequences(unique_sequences)
        self.insert_ensembl_mappings(mapping_df[ensembl_transcript_mask])
        self.insert_refseq_mappings(mapping_df[refseq_transcript_mask])
        self.insert_domains(pfam_tsv)

    def insert_sequences(self, unique_sequences):
        pfam_sequences = []
        for seq in unique_sequences:
            if m := self.PFAM_SEQUENCE_PATTERN.match(seq):
                seq = m.group(1)  # Strip off end bit
            pfam_sequences.append(PfamSequence(seq_id=seq))
        if pfam_sequences:
            print(f"Creating {len(pfam_sequences)} PFam sequences")
            PfamSequence.objects.bulk_create(pfam_sequences, batch_size=self.BULK_INSERT_SIZE, ignore_conflicts=True)

    def insert_ensembl_mappings(self, ensembl_mapping_df):
        # PFam maps Sequence <-> Ensembl Transcript ID without version
        transcripts_qs = Transcript.objects.filter(annotation_consortium=AnnotationConsortium.ENSEMBL)
        ensembl_transcript_ids = set(transcripts_qs.values_list("identifier", flat=True))

        num_missing_ensembl_transcripts = 0
        ensembl_mappings = []
        for _, row in ensembl_mapping_df.iterrows():
            transcript_id = row["identifier"]
            if transcript_id in ensembl_transcript_ids:
                seq_id = row["seq_id"]
                if m := self.PFAM_SEQUENCE_PATTERN.match(seq_id):
                    seq_id = m.group(1)  # Strip off end bit
                psi = PfamSequenceIdentifier(pfam_sequence_id=seq_id, transcript_id=transcript_id)
                ensembl_mappings.append(psi)
            else:
                num_missing_ensembl_transcripts += 1

        if ensembl_mappings:
            print(f"Inserting {len(ensembl_mappings)} ENSEMBL transcripts")
            PfamSequenceIdentifier.objects.bulk_create(ensembl_mappings, batch_size=self.BULK_INSERT_SIZE)
        if num_missing_ensembl_transcripts:
            print(f"Missing {num_missing_ensembl_transcripts} ENSEMBL transcripts")

    def insert_refseq_mappings(self, refseq_mapping_df):
        refseq_transcript_ids = set()
        refseq_transcript_version_id_by_accession = {}
        tv_qs = TranscriptVersion.objects.filter(transcript__annotation_consortium=AnnotationConsortium.REFSEQ)
        for pk, transcript_id, version in tv_qs.values_list("pk", "transcript_id", "version"):
            refseq_transcript_ids.add(transcript_id)
            accession = TranscriptVersion.get_accession(transcript_id, version)
            refseq_transcript_version_id_by_accession[accession] = pk

        num_missing_refseq_transcripts = 0
        num_missing_refseq_transcript_versions = 0

        refseq_mappings = []
        for _, row in refseq_mapping_df.iterrows():
            accession = row["identifier"]
            transcript_id, _ = TranscriptVersion.get_transcript_id_and_version(accession)
            if transcript_id not in refseq_transcript_ids:
                num_missing_refseq_transcripts += 1
                continue
            seq_id = row["seq_id"]
            if m := self.PFAM_SEQUENCE_PATTERN.match(seq_id):
                seq_id = m.group(1)  # Strip off end bit
            kwargs = {"pfam_sequence_id": seq_id,
                      "transcript_id": transcript_id}
            if transcript_version_id := refseq_transcript_version_id_by_accession.get(accession):
                kwargs["transcript_version_id"] = transcript_version_id
            else:
                num_missing_refseq_transcript_versions += 1
            refseq_mappings.append(PfamSequenceIdentifier(**kwargs))

        if refseq_mappings:
            print(f"Inserting {len(refseq_mappings)} RefSeq transcripts")
            PfamSequenceIdentifier.objects.bulk_create(refseq_mappings, batch_size=self.BULK_INSERT_SIZE)
        if num_missing_refseq_transcripts:
            print(f"Missing {num_missing_refseq_transcripts} RefSeq transcripts")
        if num_missing_refseq_transcript_versions:
            print(f"Missing {num_missing_refseq_transcript_versions} RefSeq transcript versions")

    def insert_domains(self, pfam_tsv):
        # We only want to insert the ones we have mappings for
        names = ["seq_id", "alignment_start", "alignment_end", "envelope_start", "envelope_end",
                 "hmm_acc", "hmm_name", "type", "hmm_start", "hmm_end", "hmm_length",
                 "bit_score", "e_value", "clan"]

        pfam_ids = set(Pfam.objects.all().values_list("pk", flat=True))
        mapping_df = pd.read_csv(pfam_tsv, header=None, names=names, sep='\t', skiprows=3)

        pfam_sequence_ids = set(PfamSequenceIdentifier.objects.all().values_list("pfam_sequence", flat=True))
        pfam_domains = []
        num_missing_pfam = 0
        for _, row in mapping_df.iterrows():
            seq_id = row["seq_id"]
            if seq_id not in pfam_sequence_ids:  # Not mapped to one of our genes
                continue

            hmm_acc = row["hmm_acc"]
            pfam_id = Pfam.get_pk_from_accession(hmm_acc)
            if pfam_id in pfam_ids:
                pf_domain = PfamDomains(pfam_sequence_id=seq_id,
                                        pfam_id=pfam_id,
                                        start=row["alignment_start"],
                                        end=row["alignment_end"])
                pfam_domains.append(pf_domain)
            else:
                num_missing_pfam += 1

        if pfam_domains:
            PfamDomains.objects.bulk_create(pfam_domains, batch_size=self.BULK_INSERT_SIZE)

        if num_missing_pfam:
            print(f"Missing {num_missing_pfam} Pfam entries")
