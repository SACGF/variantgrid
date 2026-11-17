import glob
import os
import re
import shlex
import subprocess

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


# #2841: stock AnnotSV's regulatory elements step holds every SV / regulatory element intersection in memory, so a
# range with many large SVs is killed for memory on every retry. Deployments run a fork with that step fixed
# (lgmgeo/AnnotSV#358) until a release carries it; its -version is upstream's, so this is how to tell them apart.
ANNOTSV_FORK_URL = "https://github.com/davmlaw/AnnotSV"
ANNOTSV_FORK_TAGS = {"3.5.8": "v3.5.8-vg1", "3.5.10": "v3.5.10-vg1"}
_ANNOTSV_FORK_MARKER = "proc recordRegulatedGenes"


def is_annotsv_fork_installed() -> bool:
    """ Whether ANNOTATION_ANNOTSV_BIN runs the fork's regulatory elements step, found where the AnnotSV launcher
        looks for its Tcl sources (<install>/share/tcl<version>/AnnotSV, else share/tcl/AnnotSV). """
    install_dir = os.path.dirname(os.path.dirname(os.path.realpath(settings.ANNOTATION_ANNOTSV_BIN)))
    for tcl_dir in sorted(glob.glob(os.path.join(install_dir, "share", "tcl*", "AnnotSV"))):
        regulatory_elements_tcl = os.path.join(tcl_dir, "AnnotSV-regulatoryelements.tcl")
        if os.path.exists(regulatory_elements_tcl):
            with open(regulatory_elements_tcl, encoding="utf-8") as f:
                return _ANNOTSV_FORK_MARKER in f.read()
    return False


def get_annotsv_fork_install_instructions() -> str:
    version = get_annotsv_command_line_version()
    tag = ANNOTSV_FORK_TAGS.get(version, f"the v<version>-vg<N> tag matching AnnotSV {version}")
    install_dir = os.path.dirname(os.path.dirname(settings.ANNOTATION_ANNOTSV_BIN))
    return (f"Install AnnotSV from {ANNOTSV_FORK_URL} ({tag}) as the owner of {install_dir}:\n"
            f"git clone --branch {tag} {ANNOTSV_FORK_URL}.git AnnotSV_src && cd AnnotSV_src\n"
            f"make PREFIX={shlex.quote(install_dir)} install\n"
            "Same -version and output as upstream, so no new pipeline version is needed. "
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

