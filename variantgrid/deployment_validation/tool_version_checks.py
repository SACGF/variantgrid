import logging
import os
import re
import subprocess
from subprocess import check_call

from django.conf import settings

from snpdb.models import SomalierConfig
from upload.vcf.vcf_preprocess import get_bcftools_tool_version


def _check_bcftools_version():
    tv = get_bcftools_tool_version(settings.BCFTOOLS_COMMAND)
    if m := re.match("^bcftools (\d+)\.(\d+).*?,", tv.version):
        major_version, minor_version = m.groups()
        return (int(major_version), int(minor_version)) >= (1, 20)
    return False


def _run_ensure_success(command_list, **kwargs):
    check_call(command_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    return True


def check_tool_versions():
    # If it returns True then everything is good, exception or anything else, bad
    tools = {
        "bcftools": _check_bcftools_version,
        "bedtools": lambda: _run_ensure_success(["bedtools", "--version"]),
        "zcat": lambda: _run_ensure_success(["zcat", "--version"]),
    }

    somalier_settings = settings.SOMALIER
    if somalier_settings["enabled"]:
        cfg = SomalierConfig()
        somalier_bin = cfg.get_annotation("command")
        tools["somalier"] = lambda: _run_ensure_success([somalier_bin])

    if madeline2_cmd := settings.PEDIGREE_MADELINE2_COMMAND:
        tools["ped_parser"] = lambda: _run_ensure_success(["ped_parser", "--version"])
        tools["madeline2"] = lambda: _run_ensure_success([madeline2_cmd, "--version"])

    if settings.LIFTOVER_BCFTOOLS_ENABLED:
        liftover_cmd = [
            "bcftools", "+liftover", "--version"
        ]
        env = os.environ.copy()
        env["BCFTOOLS_PLUGINS"] = settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR
        tools["bcftools +liftover"] = lambda: _run_ensure_success(liftover_cmd, env=env)

    tool_versions = {}
    for tool_name, func in tools.items():
        valid = False
        try:
            valid = func()
        except Exception as e:
            logging.error(e)
        tool_versions[tool_name] = valid

    return tool_versions
