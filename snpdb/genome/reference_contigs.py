"""

Need to be able to handle conversion between contig names and chromosome names

https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt

In [10]: hg_19_v
Out[10]: SequenceVariant(ac=NC_000001.10, type=g, posedit=150550916G>A)

In [11]: hg_38_v
Out[11]: SequenceVariant(ac=NC_000001.11, type=g, posedit=150578440G>A)

"""
import logging
import os

import pandas as pd
from django.conf import settings
from django.utils.text import slugify

from library.pandas_utils import read_csv_skip_header
from snpdb.models import SequenceRole, AssemblyMoleculeType

GRCH37 = "GCF_000001405.25_GRCh37.p13_assembly_report.txt"
ASSEMBLY_REPORTS = {
    "hg19": GRCH37,  # use this but then slightly modify chrM
    "GRCh37": GRCH37,
    "GRCh38": "GCF_000001405.39_GRCh38.p13_assembly_report.txt",
    "T2T-CHM13v2.0": "GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt",
}


def _get_assembly_report_df(build):
    column_names = [
        "Sequence-Name", "Sequence-Role", "Assigned-Molecule", "Assigned-Molecule-Location/Type",
        "GenBank-Accn", "Relationship", "RefSeq-Accn", "Assembly-Unit", "Sequence-Length", "UCSC-style-name"
    ]

    assembly_report = ASSEMBLY_REPORTS[build]
    assembly_filename = os.path.join(settings.BASE_DIR, "snpdb", "genome", "reference", assembly_report)
    return read_csv_skip_header(assembly_filename, sep='\t', names=column_names)


def get_assembly_report_df(build):
    df = _get_assembly_report_df(build)

    if build in ["GRCh37", "hg19"]:
        df_grch38 = _get_assembly_report_df("GRCh38")

        # The original GRCh37 didn't have MT added after genome build - but now it has it as we are using GRCh37.p13

        # hg19 is same as grch37 except for mitochondrion
        # (hg19 = NC_001807 and GRCh37 = NC_012920)
        mt_index = df_grch38["Sequence-Name"] == "MT"
        mt_grch38 = df_grch38[mt_index]

        if build == "hg19":
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
            df = pd.concat([df, mt_hg19])

    # We want MT to be at the end of the chromosomes
    # Otherwise snpEFF gives errors about contigs not being in karyotypic order
    mt_index = df["Sequence-Name"] == 'MT'
    molecules_index = df["Sequence-Role"] == 'assembled-molecule'
    other_molecules_index = molecules_index & ~mt_index

    df = pd.concat((df[other_molecules_index], df[mt_index], df[~molecules_index]))
    return df


def create_build_and_contigs(get_model: callable, build_name, alias=None, igv_genome=None):
    """ Explicitly pass in the models, because we want this to be able to work from a Django data migration """
    GenomeBuild = get_model("GenomeBuild")
    Contig = get_model("Contig")
    GenomeBuildContig = get_model("GenomeBuildContig")

    df = get_assembly_report_df(build_name)

    existing_build = GenomeBuild.objects.filter(name=build_name).first()
    if existing_build:
        raise ValueError(
            f"GenomeBuild {build_name} already exists. Use --clear if you want to delete (will delete all variants!)")

    kwargs = {"name": build_name}
    # This is so that if we create new builds, or squish migrations, we don't need the 1 off script to populate slug
    if hasattr(GenomeBuild, "slug"):
        kwargs["slug"] = slugify(build_name)

    if alias:
        kwargs["alias"] = alias
    if igv_genome:
        kwargs["igv_genome"] = igv_genome

    genome_build = GenomeBuild.objects.create(**kwargs)
    logging.info("Created build %s", genome_build)
    role_lookup = dict(((k.label, k.value) for k in SequenceRole))
    molecule_type_lookup = dict(((k.label, k.value) for k in AssemblyMoleculeType))
    molecule_type_lookup["na"] = None

    i = 0
    for _, row in df.iterrows():
        role = row["Sequence-Role"]
        molecule_type = row["Assigned-Molecule-Location/Type"]
        defaults = {
            "name": row["Sequence-Name"],
            "role": role_lookup[role],
            "assigned_molecule": row["Assigned-Molecule"],
            "molecule_type": molecule_type_lookup[molecule_type],
            "genbank_accession": row["GenBank-Accn"],
            "ucsc_name": row["UCSC-style-name"],
            "length": row["Sequence-Length"]
        }
        contig, _ = Contig.objects.get_or_create(refseq_accession=row["RefSeq-Accn"],
                                                 defaults=defaults)
        GenomeBuildContig.objects.get_or_create(genome_build=genome_build,
                                                contig=contig,
                                                defaults={"order": i})
        i += 1
