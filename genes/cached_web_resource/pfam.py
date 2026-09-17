"""
PFam entries (used for annotating protein domains). Uses data from:
 * ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
 * ftp://ftp.uniprot.org/pub/databases/uniprot/current%5Frelease/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

The domains themselves come from the InterPro API, fetched per-gene on demand - @see genes/interpro.py
"""
import ftplib
import gzip
import logging
import re
from io import BytesIO

import pandas as pd

from genes.models import (
    Pfam,
    PfamSequence,
    PfamSequenceIdentifier,
    Transcript,
    TranscriptVersion,
)
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

    # Clear existing records + recreate. This cascades away every lazily fetched PfamDomains row, which is
    # what we want - a new Pfam/UniProt release invalidates them, and they refill as genes are viewed
    Pfam.objects.all().delete()
    PfamSequence.objects.all().delete()

    num_pfam = store_pfam()
    num_sequences = store_pfam_sequences()

    cached_web_resource.description = f"{num_pfam} Pfam. {num_sequences} sequences."
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


def store_pfam_sequences() -> int:
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
    num_sequences = insert_sequences(unique_sequences)
    insert_mappings(mapping_df[ensembl_transcript_mask], AnnotationConsortium.ENSEMBL)
    insert_mappings(mapping_df[refseq_transcript_mask], AnnotationConsortium.REFSEQ)

    return num_sequences


def insert_sequences(unique_sequences) -> int:
    seq_ids = set()
    for seq in unique_sequences:
        if m := PFAM_SEQUENCE_PATTERN.match(seq):
            seq = m.group(1)  # Strip off end bit
        seq_ids.add(seq)
    pfam_sequences = [PfamSequence(seq_id=seq_id) for seq_id in seq_ids]
    if pfam_sequences:
        logging.info("Creating %d PFam sequences", len(pfam_sequences))
        PfamSequence.objects.bulk_create(pfam_sequences, batch_size=BULK_INSERT_SIZE, ignore_conflicts=True)
    return len(pfam_sequences)


def insert_mappings(mapping_df, annotation_consortium):
    ac_label = AnnotationConsortium(annotation_consortium).label
    known_transcripts = Transcript.known_transcript_ids(annotation_consortium=annotation_consortium)

    num_missing_transcripts = 0
    mappings = []
    for _, row in mapping_df.iterrows():
        accession = row["identifier"]
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(accession)
        if transcript_id not in known_transcripts:
            num_missing_transcripts += 1
            continue
        seq_id = row["seq_id"]
        if m := PFAM_SEQUENCE_PATTERN.match(seq_id):
            seq_id = m.group(1)  # Strip off end bit
        kwargs = {"pfam_sequence_id": seq_id,
                  "transcript_id": transcript_id,
                  "version": version}
        mappings.append(PfamSequenceIdentifier(**kwargs))

    if mappings:
        logging.info("Inserting %d %s transcripts", len(mappings), ac_label)
        PfamSequenceIdentifier.objects.bulk_create(mappings, batch_size=BULK_INSERT_SIZE)
    if num_missing_transcripts:
        logging.warning("Missing %d %s transcripts", num_missing_transcripts, ac_label)
