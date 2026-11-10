"""
AnnotSV's regulatory element files, and whether they are in the order AnnotSV needs (#2749).

AnnotSV intersects a run's SVs with each `*.sorted.bed` under
Annotations_Human/FtIncludedInSV/RegulatoryElements/<build>/ using `bedtools intersect -sorted`, and when that
errors falls back, without a message, to an intersect that loads the whole file (1.4 GB for the 3M-row ABC
file). The Annotations_Human 3.5 bundle ships several of them unsorted, so each deployment sorts its own copy.
The annotations directory can be read-only to VariantGrid, so nothing here writes: it finds the unsorted files
and gives the shell line an operator runs for each.

Entry points: `unsorted_regulatory_element_files` (deployment_check and the manual migration's test),
`sort_instructions` (the text both show).
"""
import glob
import os
import shlex

from django.conf import settings


def regulatory_element_files() -> list[str]:
    """ Every regulatory element file of every build this deployment runs AnnotSV on """
    filenames = []
    for annotsv_build in settings.ANNOTATION_ANNOTSV_GENOME_BUILD.values():
        build_dir = os.path.join(settings.ANNOTATION_ANNOTSV_ANNOTATIONS_DIR, "Annotations_Human", "FtIncludedInSV",
                                 "RegulatoryElements", annotsv_build)
        filenames.extend(sorted(glob.glob(os.path.join(build_dir, "*.sorted.bed"))))
    return filenames


def is_bed_sorted(filename: str) -> bool:
    """ In `LC_ALL=C sort -k1,1 -k2,2n` order, which is what `bedtools intersect -sorted` accepts alongside
        AnnotSV's own input bed: chromosomes in byte order, each appearing once, starts ascending within """
    previous = (b"", 0)
    with open(filename, "rb") as f:
        for line in f:
            chrom, start = line.split(b"\t", 2)[:2]
            key = (chrom, int(start))
            if key < previous:
                return False
            previous = key
    return True


def unsorted_regulatory_element_files() -> list[str]:
    return [filename for filename in regulatory_element_files() if not is_bed_sorted(filename)]


def sort_command(filename: str) -> str:
    """ Plain shell with absolute paths, so it runs as the files' owner on whichever machine holds them """
    quoted = shlex.quote(filename)
    quoted_tmp = shlex.quote(filename + ".tmp")
    return f"LC_ALL=C sort -k1,1 -k2,2n {quoted} > {quoted_tmp} && mv {quoted_tmp} {quoted}"


def sort_instructions(filenames: list[str]) -> str:
    commands = "\n".join(sort_command(filename) for filename in filenames)
    return (f"{len(filenames)} AnnotSV regulatory element file(s) are not sorted, so every AnnotSV run loads them "
            "into memory (about 1.4 GB more per run, #2749). AnnotSV's output is the same either way. Run these "
            "as the owner of the annotation files, on the machine that holds them:\n"
            f"{commands}\n"
            "then 'python3 manage.py deployment_check' no longer warns about them.")
