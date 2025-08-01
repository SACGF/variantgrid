import gzip
import logging
from collections import Counter
from io import TextIOWrapper

import pandas as pd
import requests
from Bio import Entrez, SeqIO

from annotation.models import CachedWebResource, GenePubMedCount
from genes.models import Gene, GeneSymbol, GeneSymbolAlias, TranscriptVersionSequenceInfoFastaFileImport, Transcript, \
    TranscriptVersionSequenceInfo, TranscriptVersion
from genes.models_enums import AnnotationConsortium, GeneSymbolAliasSource
from library.constants import MINUTE_SECS
from library.django_utils import chunked_queryset
from library.utils import sha256sum_str, iter_http_lines


def store_refseq_gene_summary_from_web(cached_web_resource: CachedWebResource):
    retrieve_refseq_gene_summaries()

    refseq_genes_qs = Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ)
    num_summaries = refseq_genes_qs.exclude(summary__isnull=True).exclude(summary='').count()
    cached_web_resource.description = f"{num_summaries} RefSeq genes w/summary"
    cached_web_resource.save()


def store_refseq_gene_info_from_web(cached_web_resource: CachedWebResource):
    GENE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz"
    r = requests.get(GENE_INFO_URL, stream=True, timeout=MINUTE_SECS)
    f = gzip.GzipFile(fileobj=r.raw)
    gene_info_df = pd.read_csv(f, sep='\t')

    known_gene_symbols = GeneSymbol.get_upper_case_lookup()
    known_symbols = set(known_gene_symbols) | set(GeneSymbolAlias.get_upper_case_lookup())

    GeneSymbolAlias.objects.filter(source=GeneSymbolAliasSource.NCBI).delete()
    symbols_synonyms = gene_info_df[["Symbol", "Synonyms"]]
    symbols_synonyms.loc[:, "Symbol"] = symbols_synonyms.loc[:, "Symbol"].str.upper()
    symbols_synonyms.loc[:, "Synonyms"] = symbols_synonyms.loc[:, "Synonyms"].str.upper()

    gene_symbols = []
    gene_symbol_aliases = []
    for _, (symbol, synonyms) in symbols_synonyms.iterrows():
        gene_symbol_id = known_gene_symbols.get(symbol)
        if gene_symbol_id is None:
            gene_symbols.append(GeneSymbol(symbol=symbol))
            gene_symbol_id = symbol

        if synonyms == "-":
            continue

        for s in synonyms.split("|"):
            if s not in known_symbols:
                s = s.upper()
                known_symbols.add(s)
                if s == gene_symbol_id:
                    continue  # Our aliases are case-insensitive so no need to store these
                gene_symbol_aliases.append(GeneSymbolAlias(alias=s,
                                                           gene_symbol_id=gene_symbol_id,
                                                           source=GeneSymbolAliasSource.NCBI))

    if gene_symbols:
        GeneSymbol.objects.bulk_create(gene_symbols, ignore_conflicts=True)
    if gene_symbol_aliases:
        GeneSymbolAlias.objects.bulk_create(gene_symbol_aliases, ignore_conflicts=True)

    cached_web_resource.description = f"{len(gene_symbol_aliases)} new gene symbol aliases"
    cached_web_resource.save()


def retrieve_refseq_gene_summaries():
    # 10k limit of return data from NCBI
    BATCH_SIZE = 2000

    refseq_genes_qs = Gene.objects.filter(annotation_consortium=AnnotationConsortium.REFSEQ, summary__isnull=True)
    refseq_genes_qs = refseq_genes_qs.exclude(pk__startswith=Gene.FAKE_GENE_ID_PREFIX)

    for genes_qs in chunked_queryset(refseq_genes_qs, BATCH_SIZE):
        gene_ids = genes_qs.values_list("pk", flat=True)
        try:
            gene_annotation = retrieve_entrez_gene_annotation(gene_ids)
        except RuntimeError:
            logging.error("API failure for gene_ids:")
            logging.error(",".join(gene_ids))
            raise

        gene_records = []
        for annot in gene_annotation:
            gene_id = annot.attributes["uid"]
            summary = annot["Summary"]
            gene_records.append(Gene(pk=gene_id, summary=summary))

        if gene_records:
            Gene.objects.bulk_update(gene_records, ["summary"])


def retrieve_entrez_gene_annotation(id_list):
    """ Annotates Entrez Gene IDs (not gene symbol) using Bio.Entrez

        Based on Jinghua (Frank) Feng's code - construct gene reference data """

    request = Entrez.epost("gene", id=",".join(id_list))
    result = Entrez.read(request)

    web_env = result["WebEnv"]
    query_key = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=web_env, query_key=query_key)
    document = Entrez.read(data)

    annotations = document["DocumentSummarySet"]["DocumentSummary"]
    return annotations


def store_refseq_sequence_info_from_web(cached_web_resource: CachedWebResource):
    SEQUENCE_INFO_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz"

    known_transcripts = set(Transcript.objects.all().values_list("identifier", flat=True))
    if not known_transcripts:
        raise ValueError("No transcripts! Insert them first!")

    r = requests.get(SEQUENCE_INFO_URL, stream=True, timeout=MINUTE_SECS)
    f = gzip.GzipFile(fileobj=r.raw)
    text_f = TextIOWrapper(f)

    annotation_consortium = AnnotationConsortium.REFSEQ

    sha256_hash = sha256sum_str(SEQUENCE_INFO_URL)
    if existing_import := TranscriptVersionSequenceInfoFastaFileImport.objects.filter(sha256_hash=sha256_hash).first():
        print(f"Deleting existing TranscriptVersionSequenceInfos for fasta import {sha256_hash}")
        existing_import.delete()

    fasta_import = TranscriptVersionSequenceInfoFastaFileImport.objects.create(sha256_hash=sha256_hash,
                                                                               annotation_consortium=annotation_consortium,
                                                                               filename=SEQUENCE_INFO_URL)
    unknown_transcripts = []
    unknown_transcript_prefixes = Counter()
    records = []
    for record in SeqIO.parse(text_f, "fasta"):
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
        TranscriptVersionSequenceInfo.objects.bulk_create(records, ignore_conflicts=True, batch_size=2000)

    cached_web_resource.description = f"{num_records} TranscriptVersionSequenceInfo records"
    cached_web_resource.save()


def store_gene2pubmed_from_web(cached_web_resource: CachedWebResource):
    print("store_gene2pubmed_from_web - starting")
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
