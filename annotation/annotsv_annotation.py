import os
import re
import subprocess

from django.conf import settings

from genes.models_enums import AnnotationConsortium
from snpdb.models.models_genome import GenomeBuild


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


def get_annotsv_tsv_filename(vcf_filename: str, output_dir: str) -> str:
    """ Path AnnotSV writes its TSV to for a given input VCF - <basename>.annotated.tsv inside output_dir.
        A pure function of the input name, so any party holding the dump path can name the TSV without
        reading it off the AnnotationRun row (#1660). """
    base = os.path.splitext(os.path.basename(vcf_filename))[0]
    return os.path.join(output_dir, f"{base}.annotated.tsv")


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
