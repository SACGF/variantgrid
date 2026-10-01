""" Checks an uploaded VCF's SNV REF bases against the build it is about to be imported as (#2030).

    The build comes from the header (contig lengths, ##reference) or upload metadata, and a header pasted from
    the wrong dictionary is otherwise trusted absolutely - bcftools norm --check-ref=s then silently rewrites
    every REF and the import looks fine while every variant sits in the wrong gene. A wrong build shows up as
    ~75% mismatches (chance), a right one as ~0, so a sample from the start of the file is enough to tell.

    Entry point: check_vcf_ref_matches_build, called by upload/vcf/vcf_preprocess.py:preprocess_vcf """

import logging
from dataclasses import dataclass
from typing import Optional

import cyvcf2
from django.conf import settings

from snpdb.models.models_genome import GenomeBuild

REF_MISMATCH_MIN_SNVS_TO_FAIL = 10  # A handful of hand-typed records with wrong REFs isn't a wrong build
SNV_BASES = frozenset("ACGT")


class VCFRefMismatchError(ValueError):
    pass


@dataclass
class RefMismatchCount:
    genome_build: GenomeBuild
    num_checked: int
    num_mismatched: int

    @property
    def fraction(self) -> float:
        return self.num_mismatched / self.num_checked if self.num_checked else 0.0

    def __str__(self):
        return f"{self.fraction:.0%} ({self.num_mismatched}/{self.num_checked}) mismatch {self.genome_build}"


def read_snvs(vcf_filename: str, max_snvs: int) -> list[tuple[str, int, str]]:
    """ (chrom, 1-based pos, REF) of the first max_snvs biallelic ACGT SNVs """
    snvs = []
    for variant in cyvcf2.Reader(vcf_filename):
        ref = variant.REF.upper()
        alts = variant.ALT
        if ref in SNV_BASES and len(alts) == 1 and alts[0].upper() in SNV_BASES:
            snvs.append((variant.CHROM, variant.POS, ref))
            if len(snvs) >= max_snvs:
                break
    return snvs


def count_ref_mismatches(snvs: list[tuple[str, int, str]], genome_build: GenomeBuild, fasta) -> RefMismatchCount:
    """ fasta is indexed fasta[chrom][start:end] (GenomeFasta.fasta). SNVs on a contig the build or fasta
        doesn't have, or on an ambiguous (non-ACGT) reference base, are not counted """
    num_checked = 0
    num_mismatched = 0
    for chrom, pos, ref in snvs:
        try:
            fasta_ref = fasta[chrom][pos - 1:pos].upper()
        except (KeyError, ValueError):
            continue
        if fasta_ref not in SNV_BASES:
            continue
        num_checked += 1
        if fasta_ref != ref:
            num_mismatched += 1
    return RefMismatchCount(genome_build, num_checked, num_mismatched)


def _count_ref_mismatches_for_build(snvs, genome_build: GenomeBuild) -> Optional[RefMismatchCount]:
    try:
        fasta = genome_build.genome_fasta.fasta
    except (FileNotFoundError, KeyError) as e:
        logging.info("VCF REF check: skipping %s - no reference fasta: %s", genome_build, e)
        return None
    return count_ref_mismatches(snvs, genome_build, fasta)


def get_ref_mismatch_message(build_count: RefMismatchCount,
                             other_build_counts: list[RefMismatchCount]) -> Optional[tuple[str, bool]]:
    """ (message, fail) when build_count is over a threshold, else None """
    fraction = build_count.fraction
    if fraction <= settings.VCF_IMPORT_REF_MISMATCH_WARN_FRACTION:
        return None

    fail = (fraction > settings.VCF_IMPORT_REF_MISMATCH_FAIL_FRACTION
            and build_count.num_checked >= REF_MISMATCH_MIN_SNVS_TO_FAIL)
    message = f"SNV REF bases: {build_count} reference."
    if other_build_counts:
        better_builds = [obc for obc in other_build_counts
                         if obc.num_checked and obc.fraction <= settings.VCF_IMPORT_REF_MISMATCH_WARN_FRACTION]
        if better_builds:
            message += " The VCF looks to have been called against " + \
                ", ".join(f"{obc.genome_build} ({obc.fraction:.0%} mismatch)" for obc in better_builds) + "."
        else:
            message += " Other builds: " + ", ".join(str(obc) for obc in other_build_counts) + "."
    if fail:
        message += " Import stopped - if the build is wrong, re-upload with 'genome_build' in the upload metadata."
    else:
        message += " Mismatched REF bases are replaced with the reference base on import."
    return message, fail


def check_vcf_ref_matches_build(vcf_filename: str, genome_build: GenomeBuild) -> Optional[str]:
    """ Raises VCFRefMismatchError when most SNV REFs disagree with genome_build's fasta, returns a warning
        message when more than a few do, else None """
    snvs = read_snvs(vcf_filename, settings.VCF_IMPORT_REF_CHECK_SNVS)
    if not snvs:
        return None

    build_count = count_ref_mismatches(snvs, genome_build, genome_build.genome_fasta.fasta)
    if build_count.fraction <= settings.VCF_IMPORT_REF_MISMATCH_WARN_FRACTION:
        return None

    other_build_counts = []
    for other_build in GenomeBuild.builds_with_annotation().exclude(pk=genome_build.pk):
        if other_count := _count_ref_mismatches_for_build(snvs, other_build):
            other_build_counts.append(other_count)

    message, fail = get_ref_mismatch_message(build_count, other_build_counts)
    if fail:
        raise VCFRefMismatchError(message)
    return message
