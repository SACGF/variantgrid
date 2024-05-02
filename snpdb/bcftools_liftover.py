import logging
import os
from subprocess import CalledProcessError

from django.conf import settings

from annotation.vep_annotation import VEPConfig
from library.utils import execute_cmd
from snpdb.models import GenomeBuild


def bcftools_liftover(source_vcf: str, source_genome_build: GenomeBuild,
                      out_vcf: str, out_genome_build: GenomeBuild) -> tuple[str, int]:
    """ returns reject file and items to process if any failed """
    vc_37 = VEPConfig(GenomeBuild.grch37())
    vc_38 = VEPConfig(GenomeBuild.grch38())


    ASSEMBLIES = {
        # These are UCSC ones
        # "GRCh37": os.path.join(settings.LIFTOVER_BCFTOOLS_DATA_DIR, "GRCh37", "human_g1k_v37.fasta"),
        # "GRCh38": os.path.join(settings.LIFTOVER_BCFTOOLS_DATA_DIR, "GRCh38", "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"),
        # Ensembl ones
        "GRCh37": vc_37["fasta"],
        "GRCh38": vc_38["fasta"],
    }

    source_fasta_filename = ASSEMBLIES.get(source_genome_build.name)
    dest_fasta_filename = ASSEMBLIES.get(out_genome_build.name)

    chain_dir = os.path.join(settings.LIFTOVER_BCFTOOLS_DATA_DIR, "chain")
    chain_files = {
        "GRCh37": {"GRCh38": os.path.join(chain_dir, "GRCh37_to_GRCh38.chain.gz")},
        "GRCh38": {"GRCh37": os.path.join(chain_dir, "GRCh38_to_GRCh37.chain.gz")},
    }
    chain_filename = chain_files.get(source_genome_build.name, {}).get(out_genome_build.name)

    required_files = {
        "source_fasta_filename": source_fasta_filename,
        "dest_fasta_filename": dest_fasta_filename,
        "chain_filename": chain_filename,
    }

    for description, filename in required_files.items():
        if filename is None:
            raise ValueError(f"BCFTools +liftover couldn't determine what '{description}' to use in liftover from {source_genome_build} -> {out_genome_build}''")
        if not os.path.exists(filename):
            raise ValueError(f"BCFTools +liftover '{description}' ('{filename}') does not exist.")

    reject_vcf = get_reject_vcf_filename(out_vcf)
    bcftools_cmd = "bcftools"  # Installed in path
    cmd = [
        bcftools_cmd,
        "+liftover",
        source_vcf,
        "--",
        "--src-fasta-ref", source_fasta_filename,
        "--fasta-ref", dest_fasta_filename,
        "--chain", chain_filename,
        "--reject", reject_vcf,
    ]

    env = {
        "BCFTOOLS_PLUGINS": settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR,
    }

    logging.info("Executing BCFTools +liftover")
    liftover_cmd_str = " ".join(cmd)
    # bcftools +liftover doesn't sort VCF - @see https://github.com/freeseek/score/issues/8
    sort_cmd = [
        bcftools_cmd,
        "sort",
        "-Oz",
        "-o", out_vcf,
    ]
    sort_cmd_str = " ".join(sort_cmd)
    cmd_str = " | ".join((liftover_cmd_str, sort_cmd_str))
    logging.info(cmd_str)
    return_code, std_out, std_err = execute_cmd([cmd_str], env=env, shell=True)
    print(f"return_code: {return_code}")
    if std_out:
        print("stdout:")
        print(std_out)
    if std_err:
        print("stderr:")
        print(std_err)

    if return_code:
        raise CalledProcessError(return_code, cmd_str, output=std_out, stderr=std_err)

    items_to_process = 0
    if os.path.exists(reject_vcf):
        items_to_process = count_non_header_lines(reject_vcf)

    return reject_vcf, items_to_process


def get_reject_vcf_filename(out_vcf_filename: str) -> str:
    basename = os.path.basename(out_vcf_filename)
    dirname = os.path.dirname(out_vcf_filename)
    return os.path.join(dirname, "rejected_" + basename)


def count_non_header_lines(filename) -> int:
    count = 0
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                count += 1
    return count
