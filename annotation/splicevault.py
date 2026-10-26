"""
SpliceVault data for RefSeq VEP runs (#2426).

The SpliceVault VEP plugin only returns a row whose transcript column equals the annotated transcript's
stable_id, and Ensembl's file is keyed by unversioned ENST. A --refseq run's stable_ids are versioned NMs
(the same ids the transcript_blocklist file holds), so the plugin matches nothing until the file is rewritten
to MANE RefSeq ids. Rows for transcripts with no MANE pair are dropped - those are all a RefSeq run loses.

Entry points: mane_ensembl_to_refseq (the id map from genes MANE), convert_splicevault_file (rewrite, bgzip,
tabix), used by `manage.py splicevault_to_refseq`.
"""
import gzip
import os
from dataclasses import dataclass

import pysam

from genes.models.models_gene import MANE

SPLICEVAULT_TRANSCRIPT_COLUMN = 4
CONVERTED_SUFFIX = ".refseq_mane.tsv.gz"


@dataclass
class SpliceVaultConversion:
    rows_read: int = 0
    rows_written: int = 0
    transcripts_converted: int = 0
    transcripts_dropped: int = 0

    @property
    def rows_dropped(self) -> int:
        return self.rows_read - self.rows_written


def mane_ensembl_to_refseq() -> dict[str, str]:
    """ Unversioned ENST (as SpliceVault keys it) -> versioned NM (as VEP's RefSeq cache names it) """
    pairs = MANE.objects.filter(ensembl_transcript_version__isnull=False, refseq_transcript_version__isnull=False)
    return {
        enst: f"{nm}.{nm_version}"
        for enst, nm, nm_version in pairs.values_list("ensembl_transcript_version__transcript_id",
                                                      "refseq_transcript_version__transcript_id",
                                                      "refseq_transcript_version__version")
    }


def default_output_filename(input_filename: str) -> str:
    basename = os.path.basename(input_filename)
    stem = basename.removesuffix(".gz").removesuffix(".tsv")
    return os.path.join(os.path.dirname(input_filename), stem + CONVERTED_SUFFIX)


def convert_splicevault_file(input_filename: str, output_filename: str,
                             ensembl_to_refseq: dict[str, str]) -> SpliceVaultConversion:
    """ Writes output_filename (bgzip) and its .tbi. Rewriting a column keeps the input's sort order. """
    conversion = SpliceVaultConversion()
    converted = set()
    dropped = set()
    with gzip.open(input_filename, "rt") as in_f, pysam.BGZFile(output_filename, "wb") as out_f:
        for line in in_f:
            if line.startswith("#"):
                out_f.write(line.encode())
                continue
            conversion.rows_read += 1
            columns = line.split("\t")
            enst = columns[SPLICEVAULT_TRANSCRIPT_COLUMN]
            if refseq := ensembl_to_refseq.get(enst):
                converted.add(enst)
                columns[SPLICEVAULT_TRANSCRIPT_COLUMN] = refseq
                out_f.write("\t".join(columns).encode())
                conversion.rows_written += 1
            else:
                dropped.add(enst)

    pysam.tabix_index(output_filename, seq_col=0, start_col=1, end_col=1, meta_char="#", force=True)
    conversion.transcripts_converted = len(converted)
    conversion.transcripts_dropped = len(dropped)
    return conversion
