import glob
import os
import pwd
import re
import shlex
import subprocess
from typing import Optional

from django.conf import settings

from genes.enums import AnnotationConsortium
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
        # Only full lines are imported (annotation/vcf_files/bulk_annotsv_tsv_inserter.py). Split lines are still
        # annotated for the ranking, but not written - nor, on our fork, held until the output is written (#2841)
        "-annotationMode", "full",
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


# #2841: stock AnnotSV holds whole regulatory element intersects, and every split line, in memory, so a range with
# many large SVs is killed for memory on every retry, and never applies ACMG criteria 4O/2C/2F from partially
# overlapping benign SVs. Deployments run a fork with those fixed until a release carries them. The fork's output
# differs from upstream's, so its -version does too (<upstream>-vg<N>, the suffix changing only with the output).
ANNOTSV_FORK_URL = "https://github.com/davmlaw/AnnotSV"
ANNOTSV_FORK_TAGS = {"3.5.8": "v3.5.8-vg4", "3.5.10": "v3.5.10-vg4"}
# Each fix in the fork: the Tcl source it changes, and a line only the fixed source has
ANNOTSV_FORK_FIXES = {
    "AnnotSV-regulatoryelements.tcl": "proc recordRegulatedGenes",  # lgmgeo/AnnotSV#358
    "AnnotSV-write.tcl": 'if {$g_AnnotSV(annotationMode) eq "full" && $AnnotationMode eq "split"} {',  # lgmgeo/AnnotSV#359
    "AnnotSV-benignsv.tcl": 'set L_benignGenes [split [join $L_genes($benign_coord) ";"] ";"]',  # lgmgeo/AnnotSV#360
}


def get_annotsv_upstream_version(version: str) -> str:
    """ The upstream release a -version is, or is a fork of: "3.5.8-vg4" -> "3.5.8" """
    return version.split("-vg")[0]


def _get_annotsv_tcl_dir() -> Optional[str]:
    """ Where the AnnotSV launcher looks for its Tcl sources: <install>/share/tcl<version>/AnnotSV, else
        share/tcl/AnnotSV """
    install_dir = os.path.dirname(os.path.dirname(os.path.realpath(settings.ANNOTATION_ANNOTSV_BIN)))
    for tcl_dir in sorted(glob.glob(os.path.join(install_dir, "share", "tcl*", "AnnotSV"))):
        if os.path.exists(os.path.join(tcl_dir, "AnnotSV-write.tcl")):
            return tcl_dir
    return None


def get_annotsv_fork_fixes_installed() -> set[str]:
    """ The ANNOTSV_FORK_FIXES sources the installed AnnotSV has the fixed version of: none for stock AnnotSV """
    installed = set()
    if tcl_dir := _get_annotsv_tcl_dir():
        for tcl_filename, marker in ANNOTSV_FORK_FIXES.items():
            tcl_path = os.path.join(tcl_dir, tcl_filename)
            if os.path.exists(tcl_path):
                with open(tcl_path, encoding="utf-8") as f:
                    if marker in f.read():
                        installed.add(tcl_filename)
    return installed


def is_annotsv_fork_installed() -> bool:
    """ Whether ANNOTATION_ANNOTSV_BIN runs the current fork, with every fix in ANNOTSV_FORK_FIXES """
    return get_annotsv_fork_fixes_installed() == set(ANNOTSV_FORK_FIXES)


def get_annotsv_fork_install_instructions() -> str:
    """ Shell to install the fork over the existing AnnotSV: code only - make install leaves the annotations bundle
        alone. A version the fork has no tag for moves to the latest tag. Either way -version changes, so the fork
        is a new pipeline version. """
    version = get_annotsv_command_line_version()
    install_dir = os.path.dirname(os.path.dirname(settings.ANNOTATION_ANNOTSV_BIN))
    tag = ANNOTSV_FORK_TAGS.get(get_annotsv_upstream_version(version))
    if tag is None:
        tag = list(ANNOTSV_FORK_TAGS.values())[-1]
        version_note = f"The fork has no tag for the installed AnnotSV {version!r}, so this moves to {tag}."
    else:
        version_note = f"AnnotSV {version} is installed, so use {tag}."
    version_note += (f" It reports -version {tag.removeprefix('v')}: afterwards register and promote it as a new "
                     "pipeline version (manage.py create_new_annotation_pipeline_version --pipeline-type=A, then "
                     "promote on the annotation runs page), which re-runs AnnotSV on every range to correct the stored ACMG classes. "
                     "AnnotSV runs are refused until it is promoted.")
    try:
        owner = pwd.getpwuid(os.stat(install_dir).st_uid).pw_name
    except (OSError, KeyError):
        owner = "<owner of the AnnotSV install>"
    backup = os.path.join(install_dir, f"code_backup_before_{tag}.tgz")
    q_install_dir = shlex.quote(install_dir)
    install_script = (f"set -e; cd {q_install_dir}; "
                      f"tar czf {shlex.quote(backup)} $(ls -d bin etc Makefile share/bash share/doc share/python3 share/tcl* 2>/dev/null); "
                      "src=$(mktemp -d); "
                      f"git clone --quiet --branch {tag} {ANNOTSV_FORK_URL}.git $src; "
                      f"make -C $src PREFIX={q_install_dir} install; rm -rf $src")
    return (f"{version_note} Only AnnotSV's code is replaced - the annotations bundle in "
            f"{settings.ANNOTATION_ANNOTSV_ANNOTATIONS_DIR} is not downloaded again. On the machine that runs AnnotSV "
            f"(the {settings.ANNOTATION_ANNOTSV_WORKER_QUEUE} celery worker), as the owner of {install_dir}:\n"
            f"sudo -u {owner} bash -c {shlex.quote(install_script)}\n"
            f"That backs up the code it replaces to {backup} (to roll back: "
            f"sudo -u {owner} tar xzf {shlex.quote(backup)} -C {q_install_dir}). 'make install' also re-clones "
            "share/python3/variantconvert and may print 'variantconvert not installed' (pip refusing the system "
            "Python) and 'Annotations need to be installed': both are harmless, VariantGrid reads AnnotSV's TSV. "
            "Check: 'python3 manage.py deployment_check' no longer warns about the AnnotSV fork. "
            "https://github.com/SACGF/variantgrid/wiki/Install-AnnotSV")


def _regulatory_elements_glob() -> str:
    return os.path.join(settings.ANNOTATION_ANNOTSV_ANNOTATIONS_DIR, "Annotations_Human", "FtIncludedInSV",
                        "RegulatoryElements", "*", "*.sorted.bed")


_SORT_CHECK = ["sort", "-c", "-s", "-k1,1", "-k2,2n"]  # -s: key order only, which is what bedtools -sorted needs


def get_annotsv_unsorted_regulatory_element_files() -> list[str]:
    """ #2749: AnnotSV intersects each of these with `bedtools intersect -sorted`, and when that errors falls
        back, without a message, to an intersect that loads the whole file (1.4 GB for the 3M-row ABC file).
        The Annotations_Human 3.5 bundle ships several of them unsorted (lgmgeo/AnnotSV#356). """
    env = {**os.environ, "LC_ALL": "C"}
    return [filename for filename in sorted(glob.glob(_regulatory_elements_glob()))
            if subprocess.run([*_SORT_CHECK, filename], env=env, capture_output=True, check=False).returncode]


def get_annotsv_sort_regulatory_elements_script() -> str:
    """ Shell that sorts the unsorted regulatory element files in place. Run by hand as the files' owner, on
        whichever machine holds them: the annotations directory is read-only to VariantGrid on some
        deployments. Only these files - re-sorting the bundle's others changes AnnotSV's output. """
    pattern = shlex.quote(os.path.dirname(os.path.dirname(_regulatory_elements_glob()))) + "/*/*.sorted.bed"
    check = " ".join(_SORT_CHECK)
    return (f"for f in {pattern}; do\n"
            f"  LC_ALL=C {check} \"$f\" 2>/dev/null || "
            "{ LC_ALL=C sort -k1,1 -k2,2n \"$f\" > \"$f.tmp\" && mv \"$f.tmp\" \"$f\" && echo \"sorted $f\"; }\n"
            "done")

