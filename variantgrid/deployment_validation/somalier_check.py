import os
import subprocess
from subprocess import check_output
from typing import Optional

from library.log_utils import log_traceback
from snpdb.models import GenomeBuild, SomalierConfig


def verify_somalier_config() -> Optional[str]:
    somalier_cfg = SomalierConfig()
    somalier_bin = somalier_cfg.get_annotation("command")
    somalier = None
    try:
        somalier_output = check_output([somalier_bin], stderr=subprocess.STDOUT)
        somalier = somalier_output.decode().split("\n", 1)[0]
    except Exception:
        log_traceback()

    return somalier


def check_somalier() -> dict:
    somalier_data = {
        "somalier_config": {
            "valid": verify_somalier_config(),
            "fix": "Install Somalier, and place in path https://github.com/brentp/somalier/"
        }
    }
    somalier_cfg = SomalierConfig()
    for genome_build in GenomeBuild.builds_with_annotation():
        sites_path = somalier_cfg.get_sites(genome_build)
        somalier_data[f"somalier_sites_file_{genome_build.name}"] = {
            "valid": os.path.exists(sites_path),
            "fix": f"Download somalier sites VCF to {sites_path} - see scripts/install/get_somalier_release.sh",
        }

        try:
            somalier_cfg.get_sites_vcf(genome_build)
            valid = True
        except Exception:
            valid = False
        sites_vcf_name = somalier_cfg.get_sites_vcf_name(genome_build)
        somalier_data[f"somalier_sites_vcf_{genome_build.name}"] = {
            "valid": valid,
            "fix": f"Upload/import '{sites_vcf_name}' as a VCF for {genome_build} (the VCF name must match the sites filename)",
        }
    return somalier_data
