"""
deployment_check "Tool versions": every external binary the server shells out to, checked the way it is
used. check_tool_versions is the entry point; check_vcf_split_pipe runs the real VCF import split stage
(GNU split --filter | bash | bgzip) on a tiny VCF, since a tool missing inside that filter only fails
at import time rather than on the command line.
"""
import gzip
import logging
import os
import re
import subprocess
import tempfile

from django.conf import settings

from pedigree.graphs.pedigree_chart import get_ped_parser_command
from snpdb.models import SomalierConfig
from upload.vcf.vcf_preprocess import get_bcftools_tool_version, get_split_vcf_command

_REQUIRED_BCFTOOLS_VERSION = (1, 20)
_INSTALL_BCFTOOLS = "https://github.com/SACGF/variantgrid/wiki/Install-bcftools-liftover"
_INSTALL_HTSLIB = "Install htslib bgzip/tabix (Debian/Ubuntu: 'apt install tabix', or from https://github.com/samtools/htslib)"
_INSTALL_COREUTILS = "Install GNU coreutils and bash (VCF import uses 'split --filter', which BSD/busybox split lacks)"
_INSTALL_PED_PARSER_MADELINE2 = "https://github.com/SACGF/variantgrid/wiki/Install-ped_parser-and-Madeline2"


def _check_bcftools_version():
    tv = get_bcftools_tool_version(settings.BCFTOOLS_COMMAND)
    if m := re.match(r"^bcftools (\d+)\.(\d+).*?,", tv.version):
        major_version, minor_version = m.groups()
        return (int(major_version), int(minor_version)) >= _REQUIRED_BCFTOOLS_VERSION
    return False


def _check_bcftools_liftover_has_write_reject():
    """Return True if `bcftools +liftover --help` mentions --write-reject."""
    env = os.environ.copy()
    env["BCFTOOLS_PLUGINS"] = settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR
    cmd = [settings.BCFTOOLS_COMMAND, "+liftover", "--help"]

    res = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True)
    write_reject_opt = "--write-reject"
    return write_reject_opt in res.stdout or write_reject_opt in res.stderr


def _run_ensure_success(command_list, **kwargs):
    subprocess.check_call(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    return True


_SPLIT_CHECK_VCF = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
_SPLIT_CHECK_RECORDS = ["1\t1\t.\tA\tC\t.\t.\t.\n", "1\t2\t.\tA\tG\t.\t.\t.\n", "1\t3\t.\tA\tT\t.\t.\t.\n"]


def check_vcf_split_pipe() -> bool:
    """ Run the VCF import split stage exactly as run_pipe does (bash, pipefail, paths via env) on a 3 record
        VCF split 2 rows per file, and read the chunks back as bgzipped VCF with the header on each """
    with tempfile.TemporaryDirectory() as tmp_dir:
        header_filename = os.path.join(tmp_dir, "header.vcf")
        records_filename = os.path.join(tmp_dir, "records.vcf")
        split_vcf_dir = os.path.join(tmp_dir, "split")
        os.mkdir(split_vcf_dir)
        with open(header_filename, "w") as f:
            f.write(_SPLIT_CHECK_VCF)
        with open(records_filename, "w") as f:
            f.writelines(_SPLIT_CHECK_RECORDS)

        split_cmd = " ".join(get_split_vcf_command("check", split_file_rows=2))
        shell_command = f"set -o pipefail; cat {records_filename} | {split_cmd}"
        env = {**os.environ, "VG_HEADER_FILE": header_filename, "VG_SPLIT_VCF_DIR": split_vcf_dir}
        subprocess.run(shell_command, shell=True, executable="/bin/bash", env=env, check=True,
                       capture_output=True, text=True)

        split_filenames = sorted(os.listdir(split_vcf_dir))
        if split_filenames != ["check00.vcf.gz", "check01.vcf.gz"]:
            raise ValueError(f"Expected 2 split files, got {split_filenames}")
        records = []
        for filename in split_filenames:
            with gzip.open(os.path.join(split_vcf_dir, filename), "rt") as f:
                lines = f.readlines()
            if "".join(lines[:2]) != _SPLIT_CHECK_VCF:
                raise ValueError(f"{filename} does not start with the VCF header: {lines[:2]}")
            records.extend(lines[2:])
        if records != _SPLIT_CHECK_RECORDS:
            raise ValueError(f"Split files did not round-trip the records: {records}")
    return True


def check_tool_versions() -> dict:
    # If it returns True then everything is good, exception or anything else, bad
    tools = {
        "bcftools": (_check_bcftools_version, _INSTALL_BCFTOOLS),
        "bedtools": (lambda: _run_ensure_success(["bedtools", "--version"]), "Install Bedtools"),
        "zcat": (lambda: _run_ensure_success([settings.BASH_ZCAT, "--version"]), "Install ZCat (gzip package)"),
        "bgzip": (lambda: _run_ensure_success(["bgzip", "--version"]), _INSTALL_HTSLIB),
        "tabix": (lambda: _run_ensure_success(["tabix", "--version"]), _INSTALL_HTSLIB),
        "split": (lambda: _run_ensure_success(["split", "--version"]), _INSTALL_COREUTILS),
        "VCF import split pipe": (check_vcf_split_pipe, f"{_INSTALL_HTSLIB}; {_INSTALL_COREUTILS}"),
    }

    somalier_settings = settings.SOMALIER
    if somalier_settings["enabled"]:
        cfg = SomalierConfig()
        somalier_bin = cfg.get_annotation("command")
        tools["somalier"] = (lambda: _run_ensure_success([somalier_bin]), "https://github.com/brentp/somalier")

    if madeline2_cmd := settings.PEDIGREE_MADELINE2_COMMAND:
        tools["ped_parser"] = (lambda: _run_ensure_success([*get_ped_parser_command(), "--version"]), _INSTALL_PED_PARSER_MADELINE2)
        tools["madeline2"] = (lambda: _run_ensure_success([madeline2_cmd, "--version"]), _INSTALL_PED_PARSER_MADELINE2)

    if settings.LIFTOVER_BCFTOOLS_ENABLED:
        liftover_cmd = [
            "bcftools", "+liftover", "--version"
        ]
        env = os.environ.copy()
        env["BCFTOOLS_PLUGINS"] = settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR
        tools["bcftools +liftover"] = (lambda: _run_ensure_success(liftover_cmd, env=env), _INSTALL_BCFTOOLS)
        tools["bcftools +liftover --write-reject"] = (_check_bcftools_liftover_has_write_reject, _INSTALL_BCFTOOLS + " - you need the dev version that has --write-reject")

    tool_versions = {}
    for tool_name, (func, instructions) in tools.items():
        valid = False
        try:
            valid = func()
        except Exception as e:
            logging.error(e)
        tool_versions[tool_name] = {
            "valid": valid,
            "fix": instructions,
        }

    return tool_versions
