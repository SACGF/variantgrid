import subprocess
from subprocess import check_output
from typing import Optional

from library.log_utils import log_traceback
from snpdb.models import SomalierConfig, GenomeBuild


def verify_somalier_config() -> Optional[str]:
    somalier_cfg = SomalierConfig()
    somalier_bin = somalier_cfg.get_annotation("command")
    somalier = None
    try:
        somalier_output = check_output([somalier_bin], stderr=subprocess.STDOUT)
        somalier = somalier_output.decode().split("\n", 1)[0]
    except:
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
        try:
            somalier_cfg.get_sites_vcf(genome_build)
            valid = True
        except Exception as e:
            valid = False
        somalier_data[f"somalier_sites_vcf_{genome_build.name}"] = {
            "valid": valid,
            "fix": somalier_cfg.get_sites_vcf_name(genome_build),
        }
    return somalier_data
