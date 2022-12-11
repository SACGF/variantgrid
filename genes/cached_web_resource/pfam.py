"""
PFam entries (used for annotating protein domains). Uses data from:
 * ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
 * ftp://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
 * ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/proteomes/9606.tsv.gz
"""
import ftplib
import gzip
import logging
import re
from io import BytesIO

import pandas as pd

from genes.models import TranscriptVersion, Pfam, Transcript, PfamSequence, PfamSequenceIdentifier, PfamDomains
from genes.models_enums import AnnotationConsortium

BULK_INSERT_SIZE = 2000
PFAM_SEQUENCE_PATTERN = re.compile(r"(.+?)-\d+$")


def store_pfam_from_web(cached_web_resource):
    """ Pfam-A.clans.tsv
      This file contains a list of all Pfam-A families that are in clans.
      The columns are: Pfam accession, clan accession, clan ID, Pfam
      ID, Pfam description. """

    if not Transcript.objects.exists():
        raise ValueError("No transcripts - you need to import them first.")

    # Clear existing records + recreate
    Pfam.objects.all().delete()
    PfamSequence.objects.all().delete()

    num_pfam = store_pfam()
    num_domains = store_pfam_sequences_and_domains()

    cached_web_resource.description = f"{num_pfam} Pfam. {num_domains} Domains."
    cached_web_resource.save()


def store_pfam() -> int:
    PFAM_CLANS_COLUMNS = ["accession", "clan_accession", "clan_id", "pfam_id", "description"]

    logging.debug("Retrieving PFam clans via FTP")
    ftp = ftplib.FTP("ftp.ebi.ac.uk")
    ftp.login("anonymous", "anonymous")
    buffer = BytesIO()
    ftp.retrbinary('RETR /pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz', buffer.write)
    buffer.seek(0)

    pfam_list = []
    with gzip.GzipFile(fileobj=buffer) as f:
        df = pd.read_csv(f, sep='\t', names=PFAM_CLANS_COLUMNS, header=None)
        for _, row in df.iterrows():
            pfam_list.append(Pfam(pk=Pfam.get_pk_from_accession(row["accession"]),
                                  pfam_id=row["pfam_id"],
                                  description=row["description"]))

    Pfam.objects.bulk_create(pfam_list)
    return len(pfam_list)


def store_pfam_sequences_and_domains() -> int:
    logging.debug("Retrieving Uniprot mappings via FTP")
    ftp = ftplib.FTP("ftp.uniprot.org")
    ftp.login("anonymous", "anonymous")
    buffer = BytesIO()
    ftp.retrbinary('RETR /pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz', buffer.write)
    buffer.seek(0)

    with gzip.GzipFile(fileobj=buffer) as f:
        mapping_df = pd.read_csv(f, header=None, names=["seq_id", "id_type", "identifier"], sep='\t')

    ensembl_transcript_mask = mapping_df["id_type"] == "Ensembl_TRS"
    refseq_transcript_mask = mapping_df["id_type"] == "RefSeq_NT"

    # Insert sequences
    unique_sequences = mapping_df[ensembl_transcript_mask | refseq_transcript_mask]["seq_id"].unique()
    insert_sequences(unique_sequences)
    insert_mappings(mapping_df[ensembl_transcript_mask], AnnotationConsortium.ENSEMBL)
    insert_mappings(mapping_df[refseq_transcript_mask], AnnotationConsortium.REFSEQ)

    logging.debug("Retrieving Pfam protein domains via FTP")
    ftp = ftplib.FTP("ftp.ebi.ac.uk")
    ftp.login("anonymous", "anonymous")
    buffer = BytesIO()
    ftp.retrbinary('RETR /pub/databases/Pfam/current_release/proteomes/9606.tsv.gz', buffer.write)
    buffer.seek(0)
    with gzip.GzipFile(fileobj=buffer) as pfam_tsv:
        return insert_domains(pfam_tsv)


def insert_sequences(unique_sequences):
    pfam_sequences = []
    for seq in unique_sequences:
        if m := PFAM_SEQUENCE_PATTERN.match(seq):
            seq = m.group(1)  # Strip off end bit
        pfam_sequences.append(PfamSequence(seq_id=seq))
    if pfam_sequences:
        print(f"Creating {len(pfam_sequences)} PFam sequences")
        PfamSequence.objects.bulk_create(pfam_sequences, batch_size=BULK_INSERT_SIZE, ignore_conflicts=True)


def insert_mappings(mapping_df, annotation_consortium):
    ac_label = AnnotationConsortium(annotation_consortium).label
    transcript_versions_by_id = TranscriptVersion.transcript_versions_by_id(annotation_consortium=annotation_consortium)

    num_missing_transcripts = 0
    num_missing_transcript_versions = 0

    mappings = []
    for _, row in mapping_df.iterrows():
        accession = row["identifier"]
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(accession)
        versions_dict = transcript_versions_by_id.get(transcript_id)
        if versions_dict is None:
            num_missing_transcripts += 1
            continue
        seq_id = row["seq_id"]
        if m := PFAM_SEQUENCE_PATTERN.match(seq_id):
            seq_id = m.group(1)  # Strip off end bit
        kwargs = {"pfam_sequence_id": seq_id,
                  "transcript_id": transcript_id}
        if transcript_version_id := versions_dict.get(version):
            kwargs["transcript_version_id"] = transcript_version_id
        else:
            num_missing_transcript_versions += 1
        mappings.append(PfamSequenceIdentifier(**kwargs))

    if mappings:
        print(f"Inserting {len(mappings)} {ac_label} transcripts")
        PfamSequenceIdentifier.objects.bulk_create(mappings, batch_size=BULK_INSERT_SIZE)
    if num_missing_transcripts:
        print(f"Missing {num_missing_transcripts} {ac_label} transcripts")
    if num_missing_transcript_versions:
        print(f"Missing {num_missing_transcript_versions} {ac_label} transcript versions")


def insert_domains(pfam_tsv) -> int:
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
        PfamDomains.objects.bulk_create(pfam_domains, batch_size=BULK_INSERT_SIZE)

    if num_missing_pfam:
        print(f"Missing {num_missing_pfam} Pfam entries")

    return len(pfam_domains)
