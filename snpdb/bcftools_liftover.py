import logging
import os
from subprocess import CalledProcessError
from typing import Optional

from django.conf import settings

from annotation.vep_annotation import VEPConfig
from library.utils import execute_cmd
from snpdb.models import GenomeBuild, VariantCoordinate


def bcftools_pre_liftover_error_check(variant_coordinate: VariantCoordinate, source_genome_build) -> Optional[str]:
    """ Return non-null (string reason) if you want to skip bcftools liftover """
    # BCFTools requires ref match
    if variant_coordinate.is_symbolic():
        if not settings.LIFTOVER_BCFTOOLS_SYMBOLIC:
            return f"Liftover of symbolic variants disabled via {settings.LIFTOVER_BCFTOOLS_SYMBOLIC=}"
    else:
        # Symbolics will pull out reference from build so always match, no point testing
        calculated_ref = variant_coordinate.calculated_reference(source_genome_build)
        if calculated_ref != variant_coordinate.ref:
            return f"Reference='{variant_coordinate.ref}' not equal to calculated ref from genome {calculated_ref}"

        if settings.LIFTOVER_BCFTOOLS_MAX_LENGTH:
            if variant_coordinate.max_sequence_length > settings.LIFTOVER_BCFTOOLS_MAX_LENGTH:
                return f"Variant max sequence length > {settings.LIFTOVER_BCFTOOLS_MAX_LENGTH=}"
    return None


def bcftools_liftover(source_vcf: str, source_genome_build: GenomeBuild,
                      out_vcf: str, out_genome_build: GenomeBuild) -> tuple[str, int]:
    """ returns reject file and items to process if any failed """

    vep_config_per_build = {}
    for genome_build in GenomeBuild.builds_with_annotation():
        vep_config_per_build[genome_build.name] = VEPConfig(genome_build)

    source_vep_config = vep_config_per_build[source_genome_build.name]
    dest_vep_config = vep_config_per_build[out_genome_build.name]

    source_fasta_filename = source_vep_config["fasta"]
    dest_fasta_filename = dest_vep_config["fasta"]

    annotation_build_config = settings.ANNOTATION[source_genome_build.name]
    chain_filename = annotation_build_config["liftover"].get(out_genome_build.name)
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
        "--write-reject",
    ]

    env = os.environ.copy()
    env["BCFTOOLS_PLUGINS"] = settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR

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
