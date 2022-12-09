import gzip
import logging
from collections import Counter

import pandas as pd
import requests

from annotation.models import CachedWebResource
from genes.models import TranscriptVersion, GeneSymbol, GeneVersion, MANE, HGNC
from genes.models_enums import MANEStatus
from snpdb.models import GenomeBuild


def store_mane_from_web(cached_web_resource: CachedWebResource):
    MANE_VERSION = "v1.0"
    MANE_URL = f"https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.{MANE_VERSION}.summary.txt.gz"
    logging.info("Retrieving %s", MANE_URL)
    r = requests.get(MANE_URL, stream=True)
    f = gzip.GzipFile(fileobj=r.raw)

    MANE.objects.all().delete()
    genome_build = GenomeBuild.grch38()  # MANE is for GRCh38 only
    ncbi_gene_prefix = "GeneID:"
    ncbi_gene_prefix_length = len(ncbi_gene_prefix)
    uc_known_gene_symbols = GeneSymbol.get_upper_case_lookup()
    hgnc_ids_by_accession = HGNC.id_by_accession()
    gv_ids_by_accession = GeneVersion.id_by_accession(genome_build=genome_build)
    tv_ids_by_accession = TranscriptVersion.id_by_accession(genome_build=genome_build)
    mane_status_lookup = {v: k for k, v in MANEStatus.choices}

    if not all([uc_known_gene_symbols, hgnc_ids_by_accession, gv_ids_by_accession, tv_ids_by_accession]):
        raise ValueError("Need to insert genes/transcripts and HGNC")

    df = pd.read_csv(f, sep='\t')
    new_symbols = []
    records = []
    unmatched = Counter()
    for _, row in df.iterrows():
        ncbi_gene_id = row["#NCBI_GeneID"][ncbi_gene_prefix_length:]
        symbol = row["symbol"]
        if symbol.upper() not in uc_known_gene_symbols:
            new_symbols.append(GeneSymbol(pk=symbol))
        kwargs = {
            "ncbi_gene_version_id": gv_ids_by_accession.get(ncbi_gene_id),
            "ensembl_gene_version_id": gv_ids_by_accession.get(row["Ensembl_Gene"]),
            "hgnc_id": hgnc_ids_by_accession.get(row["HGNC_ID"]),
            "symbol_id": symbol,
            "refseq_transcript_version_id": tv_ids_by_accession.get(row["RefSeq_nuc"]),
            "ensembl_transcript_version_id": tv_ids_by_accession.get(row["Ensembl_nuc"]),
            "status": mane_status_lookup[row["MANE_status"]]
        }
        for k, v in kwargs.items():
            if v is None:
                unmatched[k] += 1
        records.append(MANE(**kwargs))

    if new_symbols:
        logging.info("Creating %d new gene symbols", len(new_symbols))
        GeneSymbol.objects.bulk_create(new_symbols, batch_size=2000)

    if records:
        logging.info("Inserting %d MANE records", len(records))
        MANE.objects.bulk_create(records, batch_size=2000)

    logging.info("Unmatched records: %s", unmatched)

    cached_web_resource.description = f"{MANE_VERSION}: {len(records)} records"
    cached_web_resource.save()
