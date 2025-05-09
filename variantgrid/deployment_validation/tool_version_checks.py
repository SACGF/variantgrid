import logging
import os
import re
import subprocess

from django.conf import settings

from snpdb.models import SomalierConfig
from upload.vcf.vcf_preprocess import get_bcftools_tool_version

_REQUIRED_BCFTOOLS_VERSION = (1, 20)
_INSTALL_BCFTOOLS = "https://github.com/SACGF/variantgrid/wiki/Install-bcftools-liftover"
_INSTALL_PED_PARSER_MADELINE2 = "https://github.com/SACGF/variantgrid/wiki/Install-ped_parser-and-Madeline2"


def _check_bcftools_version():
    tv = get_bcftools_tool_version(settings.BCFTOOLS_COMMAND)
    if m := re.match(r"^bcftools (\d+)\.(\d+).*?,", tv.version):
        major_version, minor_version = m.groups()
        return (int(major_version), int(minor_version)) >= _REQUIRED_BCFTOOLS_VERSION
    return False


def _check_bcftools_liftover_has_write_reject():
    env = os.environ.copy()
    env["BCFTOOLS_PLUGINS"] = settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR

    p = subprocess.Popen([settings.BCFTOOLS_COMMAND, "+liftover", "--help"], env=env,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    stdout = stdout.decode()
    stderr = stderr.decode()
    if p.returncode:
        raise subprocess.CalledProcessError(
            f"Error running '{settings.BCFTOOLS_COMMAND} --version': {stderr=}, returned: {p.returncode}")

    # They may change whether help is in stderr or stdout
    write_reject_opt = "--write-reject"
    return write_reject_opt in stdout or write_reject_opt in stderr


def _run_ensure_success(command_list, **kwargs):
    subprocess.check_call(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    return True


def check_tool_versions() -> dict:
    # If it returns True then everything is good, exception or anything else, bad
    tools = {
        "bcftools": (_check_bcftools_version, _INSTALL_BCFTOOLS),
        "bedtools": (lambda: _run_ensure_success(["bedtools", "--version"]), "Install Bedtools"),
        "zcat": (lambda: _run_ensure_success(["zcat", "--version"]), "Install ZCat"),
    }

    somalier_settings = settings.SOMALIER
    if somalier_settings["enabled"]:
        cfg = SomalierConfig()
        somalier_bin = cfg.get_annotation("command")
        tools["somalier"] = (lambda: _run_ensure_success([somalier_bin]), "https://github.com/brentp/somalier")

    if madeline2_cmd := settings.PEDIGREE_MADELINE2_COMMAND:
        tools["ped_parser"] = (lambda: _run_ensure_success(["ped_parser", "--version"]), _INSTALL_PED_PARSER_MADELINE2)
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
