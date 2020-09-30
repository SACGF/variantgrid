import logging
import os
from subprocess import CalledProcessError

from django.conf import settings

from library.utils import execute_cmd
from snpdb.models import GenomeBuild


def ncbi_remap(source_vcf: str, source_genome_build: GenomeBuild, out_vcf: str, out_genome_build: GenomeBuild):
    REMAP_SCRIPT = os.path.join(settings.BASE_DIR, "scripts", "remap_api.pl")
    report_filename = os.path.join(os.path.dirname(out_vcf), "ncbi_remap.tsv")

    # Using patch versions, I got 3% of GRCh38 -> GRCh37 liftover going to patches on first pass.
    # Using primary assembly we avoid the patches/alt loci etc.
    ASSEMBLIES = {
        "GRCh37": "GCF_000001405.13",
        "GRCh38": "GCF_000001405.26",
    }
    source_accession = ASSEMBLIES.get(source_genome_build.name, source_genome_build.accession)
    dest_accession = ASSEMBLIES.get(out_genome_build.name, out_genome_build.accession)

    cmd = [
        "perl", REMAP_SCRIPT,
        "--mode", "asm-asm",
        "--from", source_accession,
        "--dest", dest_accession,
        "--annotation", source_vcf,
        "--annot_out", out_vcf,
        "--allowdupes", "off",
        "--report_out", report_filename
    ]

    if settings.LIFTOVER_NCBI_REMAP_PERLBREW_RUNNER_SCRIPT:
        cmd.insert(0, settings.LIFTOVER_NCBI_REMAP_PERLBREW_RUNNER_SCRIPT)

    logging.info("Executing NCBI remap:")
    cmd_str = " ".join(cmd)
    logging.info(cmd_str)
    return_code, std_out, std_err = execute_cmd(cmd)
    print(f"return_code: {return_code}")
    if std_out:
        print("stdout:")
        print(std_out)
    if std_err:
        print("stderr:")
        print(std_err)

    if return_code:
        raise CalledProcessError(return_code, cmd_str, output=std_out, stderr=std_err)
