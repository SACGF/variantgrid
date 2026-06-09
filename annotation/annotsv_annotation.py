import logging
import os
import re
import subprocess

from django.conf import settings

from genes.models_enums import AnnotationConsortium
from snpdb.models.models_genome import GenomeBuild


class AnnotSVVersionMismatchError(ValueError):
    pass


def get_annotsv_command(vcf_filename: str, output_dir: str,
                        genome_build: GenomeBuild,
                        annotation_consortium: str) -> list[str]:
    build_arg = settings.ANNOTATION_ANNOTSV_GENOME_BUILD[genome_build.name]
    tx_arg = "RefSeq" if annotation_consortium == AnnotationConsortium.REFSEQ else "ENSEMBL"
    cmd = [
        settings.ANNOTATION_ANNOTSV_BIN,
        "-SVinputFile", vcf_filename,
        "-outputDir", output_dir,
        "-genomeBuild", build_arg,
        "-annotationsDir", settings.ANNOTATION_ANNOTSV_ANNOTATIONS_DIR,
        "-tx", tx_arg,
        "-SVinputInfo", "1",
        "-includeCI", "0",
        "-overwrite", "1",
    ]
    cmd.extend(settings.ANNOTATION_ANNOTSV_EXTRA_ARGS)
    return cmd


def run_annotsv(vcf_filename: str, output_dir: str,
                genome_build: GenomeBuild,
                annotation_consortium: str) -> tuple[str, int, str, str]:
    """ Returns (tsv_filename, return_code, stdout, stderr).

        Best-effort: caller logs failure on AnnotationRun. AnnotSV writes
        <basename>.annotated.tsv inside output_dir. """
    os.makedirs(output_dir, exist_ok=True)
    cmd = get_annotsv_command(vcf_filename, output_dir, genome_build, annotation_consortium)
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True, text=True,
            timeout=settings.ANNOTATION_ANNOTSV_TIMEOUT_SECONDS,
            check=False,
        )
        rc = proc.returncode
        stdout = proc.stdout or ""
        stderr = proc.stderr or ""
    except subprocess.TimeoutExpired as e:
        rc = -1
        stdout = (e.stdout or "") if isinstance(e.stdout, str) else ""
        stderr = f"AnnotSV timed out after {settings.ANNOTATION_ANNOTSV_TIMEOUT_SECONDS}s"
        logging.warning(stderr)

    base = os.path.splitext(os.path.basename(vcf_filename))[0]
    tsv_filename = os.path.join(output_dir, f"{base}.annotated.tsv")
    if rc != 0 or not os.path.exists(tsv_filename):
        logging.warning("AnnotSV failed (rc=%s): %s", rc, stderr[:2000])
    return tsv_filename, rc, stdout, stderr


def get_annotsv_command_line_version() -> str:
    """ Run `AnnotSV -version` and return the version string (eg "3.5.8"). """
    proc = subprocess.run(
        [settings.ANNOTATION_ANNOTSV_BIN, "-version"],
        capture_output=True, text=True, check=False,
    )
    out = (proc.stdout or "") + (proc.stderr or "")
    if m := re.search(r"AnnotSV\s+([\w.\-]+)", out):
        return m.group(1)
    return out.strip()


def annotsv_check_command_line_version_match(variant_annotation_version):
    """ Counterpart to vep_check_command_line_version_match. Skipped entirely
        when ANNOTATION_ANNOTSV_ENABLED is False. """
    if not settings.ANNOTATION_ANNOTSV_ENABLED:
        return

    expected_code = variant_annotation_version.annotsv_code
    if expected_code:
        actual_code = get_annotsv_command_line_version()
        if expected_code != actual_code:
            raise AnnotSVVersionMismatchError(
                f"AnnotSV out of sync! VariantAnnotationVersion has annotsv_code='{expected_code}', "
                f"installed binary reports '{actual_code}'"
            )

    expected_bundle = variant_annotation_version.annotsv_bundle
    actual_bundle = settings.ANNOTATION_ANNOTSV_BUNDLE_VERSION
    if expected_bundle and actual_bundle and expected_bundle != actual_bundle:
        raise AnnotSVVersionMismatchError(
            f"AnnotSV bundle out of sync! VariantAnnotationVersion has annotsv_bundle='{expected_bundle}', "
            f"settings.ANNOTATION_ANNOTSV_BUNDLE_VERSION='{actual_bundle}'"
        )
