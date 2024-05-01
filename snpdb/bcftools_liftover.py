import logging
import os
from subprocess import CalledProcessError

from django.conf import settings

from annotation.vep_annotation import VEPConfig
from library.utils import execute_cmd
from snpdb.models import GenomeBuild


def bcftools_liftover(source_vcf: str, source_genome_build: GenomeBuild, out_vcf: str, out_genome_build: GenomeBuild):
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

    bcftools_cmd = "bcftools"  # Installed in path

    cmd = [
        bcftools_cmd,
        "+liftover",
        "-Oz",
        "-o", out_vcf,
        source_vcf,
        "--",
        "--src-fasta-ref", source_fasta_filename,
        "--fasta-ref", dest_fasta_filename,
        "--chain", chain_filename,
    ]

    env = {
        "BCFTOOLS_PLUGINS": settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR,
    }

    logging.info("Executing BCFTools +liftover")
    cmd_str = " ".join(cmd)
    logging.info(cmd_str)
    return_code, std_out, std_err = execute_cmd(cmd, env=env)
    print(f"return_code: {return_code}")
    if std_out:
        print("stdout:")
        print(std_out)
    if std_err:
        print("stderr:")
        print(std_err)

    if return_code:
        raise CalledProcessError(return_code, cmd_str, output=std_out, stderr=std_err)
