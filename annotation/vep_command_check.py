"""
The VEP command line each VariantAnnotationVersion annotates with, per pipeline type, as a check that the
version's runs stay alike. The pins on the version (annotation/vep_config.py:vep_component_version_kwargs) halt
annotation when a setting or data file changes; this also catches what they can't see - a code change to
get_vep_command, such as a plugin's options - and reports it without halting, since most such changes are
deliberate. Entry points: check_vep_command (each VEP run) and vep_command_diffs (mismatch reports).
"""
import difflib
import logging
import os
import shlex

from django.conf import settings
from django.db import transaction

from annotation.enums import VariantAnnotationPipelineType
from annotation.models.models_version import VariantAnnotationVersion
from annotation.vep_annotation import get_vep_command
from eventlog.models import create_event
from library.enums.log_level import LogLevel
from library.log_utils import report_message

# Arguments (each followed by its value) that differ between runs or boxes without changing VEP's output
_RUN_SPECIFIC_ARGS = {"-i", "-o", "--skipped_variants_file", "--fork", "--buffer_size", "--dir_cache",
                      "--dir_plugins", "--compress_output"}

# (version pk, pipeline type, command) already reported by this process
_reported_changes: set[tuple[int, str, str]] = set()


def normalise_vep_command(cmd: list[str]) -> str:
    """ The command with what is particular to one run or box removed: wrappers ahead of VEP (prlimit,
        perlbrew), input/output files, parallelism, and the base directories data paths sit under """
    vep_binary = os.path.join(settings.ANNOTATION_VEP_CODE_DIR, "vep")
    args = cmd[cmd.index(vep_binary) + 1:]
    kept = ["vep"]
    skip_value = False
    for arg in args:
        if skip_value:
            skip_value = False
        elif arg in _RUN_SPECIFIC_ARGS:
            skip_value = True
        else:
            for base_dir in (settings.ANNOTATION_VEP_BASE_DIR, settings.ANNOTATION_BASE_DIR):
                arg = arg.replace(base_dir.rstrip("/") + "/", "")
            kept.append(arg)
    return shlex.join(kept)


def _diff(old_command: str, new_command: str) -> str:
    old_args = shlex.split(old_command)
    new_args = shlex.split(new_command)
    return "\n".join(line for line in difflib.unified_diff(old_args, new_args, "recorded", "now", lineterm="", n=1))


def check_vep_command(variant_annotation_version: VariantAnnotationVersion,
                      pipeline_type: VariantAnnotationPipelineType, cmd: list[str]):
    """ Record the command the first time a version runs this pipeline type, and report (once per process)
        when a later run's differs """
    command = normalise_vep_command(cmd)
    with transaction.atomic():
        vav = VariantAnnotationVersion.objects.select_for_update().only("vep_commands").get(
            pk=variant_annotation_version.pk)
        vep_commands = vav.vep_commands or {}
        recorded = vep_commands.get(pipeline_type)
        if recorded is None:
            vep_commands[pipeline_type] = command
            VariantAnnotationVersion.objects.filter(pk=vav.pk).update(vep_commands=vep_commands)
    variant_annotation_version.vep_commands = vep_commands

    if recorded is None or recorded == command:
        return
    report_key = (variant_annotation_version.pk, pipeline_type, command)
    if report_key in _reported_changes:
        return
    _reported_changes.add(report_key)
    pipeline_label = VariantAnnotationPipelineType(pipeline_type).label
    details = f"{variant_annotation_version} {pipeline_label}:\n{_diff(recorded, command)}"
    logging.warning("VEP command line changed: %s", details)
    report_message(message=f"VEP command line changed for {variant_annotation_version}", level="warning",
                   extra_data={"output": details})
    create_event(None, "VEP command line changed", details, severity=LogLevel.WARNING)


def vep_command_diffs(variant_annotation_version: VariantAnnotationVersion) -> list[str]:
    """ How the command each recorded pipeline type would run now differs from what the version ran """
    diffs = []
    for pipeline_type, recorded in (variant_annotation_version.vep_commands or {}).items():
        cmd = get_vep_command("in.vcf", "out.vcf", variant_annotation_version.genome_build,
                              variant_annotation_version.annotation_consortium, pipeline_type,
                              variant_annotation_version=variant_annotation_version)
        command = normalise_vep_command(cmd)
        if command != recorded:
            pipeline_label = VariantAnnotationPipelineType(pipeline_type).label
            diffs.append(f"{pipeline_label} VEP command:\n{_diff(recorded, command)}")
    return diffs
