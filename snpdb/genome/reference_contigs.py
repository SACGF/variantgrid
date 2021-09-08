"""

Need to be able to handle conversion between contig names and chromosome names

https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt

In [10]: hg_19_v
Out[10]: SequenceVariant(ac=NC_000001.10, type=g, posedit=150550916G>A)

In [11]: hg_38_v
Out[11]: SequenceVariant(ac=NC_000001.11, type=g, posedit=150578440G>A)

"""
import os

import pandas as pd
from django.conf import settings

from library.pandas_utils import read_csv_skip_header

GRCH37 = "GCF_000001405.25_GRCh37.p13_assembly_report.txt"
ASSEMBLY_REPORTS = {"hg19": GRCH37,  # alias
                    "GRCh37": GRCH37,
                    "GRCh38": "GCF_000001405.39_GRCh38.p13_assembly_report.txt"}


def _get_assembly_report_df(build):
    column_names = ["Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type", "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name"]

    assembly_report = ASSEMBLY_REPORTS[build]
    assembly_filename = os.path.join(settings.BASE_DIR, "snpdb", "genome", "reference", assembly_report)
    return read_csv_skip_header(assembly_filename, sep='\t', names=column_names)


def get_assembly_report_df(build):
    df = _get_assembly_report_df(build)

    if build in ["GRCh37", "hg19"]:
        df_grch38 = _get_assembly_report_df("GRCh38")

        # GRCh37 assembly report doesn't have MT it's added after genome build
        # hg19 is same as grch37 except for mitochondrion
        # (hg19 = NC_001807 and GRCh37 = NC_012920)
        mt_index = df_grch38["Sequence-Name"] == "MT"
        mt_grch38 = df_grch38[mt_index]

        if build == "GRCh37":
            # GRCh37:
            # Our b37 fasta.fai MT has length 16569 which is the same as NC_012920 in hg38
            # so just use that one. See https://www.ncbi.nlm.nih.gov/nuccore/251831106

            mt_grch37 = mt_grch38.copy()
            mt_grch37["UCSC-style-name"] = None  # GRCh37<->hg19 MT are different
            df = df.append(mt_grch37)

        elif build == "hg19":
            # See http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/
            # Note on chrM:
            # Since the release of the UCSC hg19 assembly, the Homo sapiens mitochondrion
            # sequence (represented as "chrM" in the Genome Browser) has been replaced in
            # GenBank with the record NC_012920. We have not replaced the original sequence,
            # NC_001807, in the hg19 Genome Browser. We plan to use the Revised Cambridge
            # Reference Sequence (rCRS, http://mitomap.org/bin/view.pl/MITOMAP/HumanMitoSeq)
            # in the next human assembly release.

            mt_hg19 = mt_grch38.copy()
            mt_hg19["RefSeq-Accn"] = "NC_001807"
            mt_hg19["GenBank-Accn"] = None
            # Length chrM 16571 from
            # http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
            mt_hg19["Sequence-Length"] = 16571
            df = df.append(mt_hg19)

    # We want MT to be at the end of the chromosomes
    # Otherwise snpEFF gives errors about contigs not being in karyotypic order
    mt_index = df["Sequence-Name"] == 'MT'
    molecules_index = df["Sequence-Role"] == 'assembled-molecule'
    other_molecules_index = molecules_index & ~mt_index

    df = pd.concat((df[other_molecules_index], df[mt_index], df[~molecules_index]))
    return df
