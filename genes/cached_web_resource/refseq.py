import logging
from collections import Counter

from Bio.SeqIO.FastaIO import SimpleFastaParser
from django.db import transaction

from annotation.models import GenePubMedCount
from genes.enums import AnnotationConsortium
from genes.models import (
    Gene,
    Transcript,
    TranscriptVersion,
    TranscriptVersionSequenceInfo,
    TranscriptVersionSequenceInfoFastaFileImport,
)
from library.constants import MINUTE_SECS
from library.utils import sha256sum_str
from library.utils.misc_utils import iter_http_lines
from snpdb.models.models_cached_web_resource import CachedWebResource

REFSEQ_SEQUENCE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz"
# The RefSeq "latest" file above only has the current version of each transcript. Historical classifications
# often cite older versions, which would otherwise fall back to slow per-transcript NCBI Entrez fetches during
# variant re-matching. We host a supplementary FASTA of those older versions and load it straight after.
REFSEQ_SUPPLEMENTARY_SEQUENCE_INFO_URL = "https://variantgrid.com/download/annotation/fasta/refseq_transcript_sequence_supplementary.fasta.gz"


def store_transcript_sequence_info_fasta(text_lines, filename: str, annotation_consortium) -> int:
    """ Parse a FASTA of transcript sequences (>accession, then sequence) and bulk-load
        TranscriptVersionSequenceInfo. Re-loading the same filename replaces its prior import, once
        the whole file has been read. text_lines is a text file handle or any iterable of lines """
    known_transcripts = set(Transcript.objects.values_list("identifier", flat=True))

    fasta_import = TranscriptVersionSequenceInfoFastaFileImport(sha256_hash=sha256sum_str(filename),
                                                                annotation_consortium=annotation_consortium,
                                                                filename=filename)
    unknown_transcripts = []
    unknown_transcript_prefixes = Counter()
    records = []
    non_blank_lines = (line for line in text_lines if line.strip())
    for title, sequence in SimpleFastaParser(non_blank_lines):
        accession = title.split(maxsplit=1)[0]
        transcript_id, version = TranscriptVersion.get_transcript_id_and_version(accession)
        if version is None:
            logging.warning("Skipping versionless accession %s", accession)
            continue
        if transcript_id not in known_transcripts:
            if transcript_id.startswith("X"):
                continue  # We don't want these
            known_transcripts.add(transcript_id)
            unknown_transcript_prefixes[transcript_id.split("_")[0]] += 1
            unknown_transcripts.append(Transcript(identifier=transcript_id,
                                                  annotation_consortium=annotation_consortium))

        records.append(TranscriptVersionSequenceInfo(transcript_id=transcript_id, version=version,
                                                     fasta_import=fasta_import,
                                                     sequence=sequence, length=len(sequence)))

    with transaction.atomic():
        existing_qs = TranscriptVersionSequenceInfoFastaFileImport.objects.filter(sha256_hash=fasta_import.sha256_hash)
        if existing_qs.exists():
            logging.info("Deleting existing TranscriptVersionSequenceInfos for fasta import %s", filename)
            existing_qs.delete()
        fasta_import.save()

        if unknown_transcripts:
            logging.info("Inserting %d unknown_transcripts (%s)", len(unknown_transcripts), unknown_transcript_prefixes)
            Transcript.objects.bulk_create(unknown_transcripts, batch_size=2000, ignore_conflicts=True)

        num_records = len(records)
        if num_records:
            TranscriptVersionSequenceInfo.objects.bulk_create(records, ignore_conflicts=True, batch_size=2000)
    return num_records


def _load_transcript_sequence_info_from_url(url: str, annotation_consortium) -> int:
    return store_transcript_sequence_info_fasta(iter_http_lines(url, timeout=MINUTE_SECS), url, annotation_consortium)


def store_refseq_sequence_info_from_web(cached_web_resource: CachedWebResource):
    known_transcripts = set(Transcript.objects.values_list("identifier", flat=True))
    if not known_transcripts:
        raise ValueError("No transcripts! Insert them first!")

    annotation_consortium = AnnotationConsortium.REFSEQ
    num_records = _load_transcript_sequence_info_from_url(REFSEQ_SEQUENCE_INFO_URL, annotation_consortium)

    # Supplementary (older-version) sequences are an optimisation - don't fail the whole load if unavailable
    supplementary = 0
    try:
        supplementary = _load_transcript_sequence_info_from_url(REFSEQ_SUPPLEMENTARY_SEQUENCE_INFO_URL,
                                                                annotation_consortium)
    except Exception as e:
        logging.warning("Could not load supplementary transcript sequences from %s: %s",
                        REFSEQ_SUPPLEMENTARY_SEQUENCE_INFO_URL, e)

    cached_web_resource.description = f"{num_records} + {supplementary} supplementary TranscriptVersionSequenceInfo records"
    cached_web_resource.save()


def store_gene2pubmed_from_web(cached_web_resource: CachedWebResource):
    logging.info("store_gene2pubmed_from_web - starting")
    homo_sapiens_tax_id = "9606"
    expected_cols = ["#tax_id", "GeneID", "PubMed_ID"]
    url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
    found_human = False
    found_header = False
    gene_counts = Counter()

    known_gene_ids = Gene.known_gene_ids(annotation_consortium=AnnotationConsortium.REFSEQ)
    if not known_gene_ids:
        raise ValueError("No RefSeq GeneIDs found.")

    skipped_genes = set()
    for line in iter_http_lines(url, timeout=60):
        if not found_header:
            # Make sure header is correct
            header_cols = line.split("\t")
            if header_cols != expected_cols:
                raise ValueError(f"Expected columns '{expected_cols}' but got '{header_cols}'")
            found_header = True

        tax_id, gene_id, _pubmed_id = line.strip().split("\t")
        if tax_id != homo_sapiens_tax_id:
            if found_human:
                break  # Sorted so no more human
            continue

        if gene_id not in known_gene_ids:
            skipped_genes.add(gene_id)
            continue

        found_human = True
        gene_counts[gene_id] += 1

    if gene_counts:
        gene_pubmed_count_records = []
        for gene_id, count in gene_counts.items():
            gpmc = GenePubMedCount(gene_id=gene_id, count=count, cached_web_resource=cached_web_resource)
            gene_pubmed_count_records.append(gpmc)
        GenePubMedCount.objects.bulk_create(gene_pubmed_count_records, batch_size=2000)

    cached_web_resource.description = f"GenePubMedCount {len(gene_counts)} genes, total count={gene_counts.total()}"
    cached_web_resource.save()

    if skipped_genes:
        examples = ", ".join(list(skipped_genes)[:10])
        logging.info("Skipped %d genes: example: %s", len(skipped_genes), examples)

    logging.info("store_gene2pubmed_from_web - Finished: %s", cached_web_resource.description)
