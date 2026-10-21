import csv
import gzip
import io
import os
import subprocess
import tempfile
from collections import Counter
from subprocess import check_output
from typing import Optional

from bgzip import BGZipWriter
from django.conf import settings

from library.genomics.vcf_writer import VCFWriter
from library.log_utils import log_traceback
from library.utils import execute_cmd
from snpdb.models import GenomeBuild, SomalierConfig
from snpdb.variants_to_vcf import (
    SOMALIER_AD_FORMAT,
    SOMALIER_ALLELE_ORDER_NOTE,
    SOMALIER_GT_FORMAT,
)
from snpdb.vcf_export_utils import get_vcf_header_from_contigs

# Enough of each allele order that a miscount can't be a coincidence, few enough that extract is instant
ALLELE_ORDER_SITES_PER_ORDER = 10
ALLELE_ORDER_DEPTH = 80
ALLELE_ORDER_SAMPLE = "allele_order_check"
ALLELE_ORDER_UNKNOWN_SAMPLE = "allele_order_unknown"
SOMALIER_ALLELE_ORDER_ISSUE = "https://github.com/brentp/somalier/issues/163"
RELATE_SITES_OPTION = "--sites"
# What a pre-0.3.5 relate says when we hand it the sites VCF the fixed one needs
RELATE_NO_SITES_OPTION = f"unknown option: {RELATE_SITES_OPTION}"


def verify_somalier_config() -> Optional[str]:
    somalier_cfg = SomalierConfig()
    somalier_bin = somalier_cfg.get_annotation("command")
    somalier = None
    try:
        somalier_output = check_output([somalier_bin], stderr=subprocess.STDOUT)
        somalier = somalier_output.decode().split("\n", 1)[0]
    except Exception:
        log_traceback()

    return somalier


def _write_allele_order_vcf(genome_build: GenomeBuild, sites_filename: str, vcf_filename: str,
                            compensate: bool) -> int:
    """ One sample hom-alt at every site and one with no call, half of the sites ones whose ALT sorts
        first - written the way snpdb.variants_to_vcf writes a somalier export with or without the
        compensation. Returns the number of records. """
    sites = []
    per_order = Counter()
    with gzip.open(sites_filename, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, position, _vcf_id, ref, alt = line.split("\t")[:5]
            # bucket on the site's own allele order, not on what we do about it, or turning the
            # compensation off would leave us looking for sites we never pick
            alt_sorts_first = alt < ref
            if per_order[alt_sorts_first] < ALLELE_ORDER_SITES_PER_ORDER:
                per_order[alt_sorts_first] += 1
                sites.append((chrom, int(position), ref, alt))  # keep the file's coordinate order
            if all(per_order[o] == ALLELE_ORDER_SITES_PER_ORDER for o in (False, True)):
                break

    header_lines = get_vcf_header_from_contigs(genome_build,
                                               samples=[ALLELE_ORDER_SAMPLE, ALLELE_ORDER_UNKNOWN_SAMPLE],
                                               use_accession=False,
                                               formats=[SOMALIER_GT_FORMAT, SOMALIER_AD_FORMAT],
                                               top_lines=[SOMALIER_ALLELE_ORDER_NOTE])
    with open(vcf_filename, "wb") as raw:
        with BGZipWriter(raw) as bgzip_f:
            f = io.TextIOWrapper(bgzip_f, encoding="utf-8", write_through=True)
            writer = VCFWriter(f, header_lines)
            for chrom, position, ref, alt in sites:
                # every call is hom-alt, so all the depth sits on whichever allele we put second
                flipped = compensate and alt < ref
                allele_depths = f"{ALLELE_ORDER_DEPTH},0" if flipped else f"0,{ALLELE_ORDER_DEPTH}"
                writer.write_record(chrom, position, ref, alt, fmt="GT:AD",
                                    sample_calls=[f"1/1:{allele_depths}", "./.:."])
            f.flush()
            f.detach()

    return len(sites)


def _somalier_ibs0(cfg: SomalierConfig, genome_build: GenomeBuild, work_dir: str, vcf_filename: str,
                   compensate: bool) -> int:
    """ extract + relate the pair with --unknown, as the all-samples relate runs, and hand back how many
        sites somalier saw them as opposite homozygotes """
    somalier_bin = cfg.get_annotation("command")
    relate_sites_args = [] if compensate else [RELATE_SITES_OPTION, cfg.get_sites(genome_build)]
    extract_dir = os.path.join(work_dir, "extract")
    relate_dir = os.path.join(work_dir, "relate")
    for d in (extract_dir, relate_dir):
        os.makedirs(d)

    extract_cmd = [somalier_bin, "extract", "--out-dir", extract_dir,
                   "--sites", cfg.get_sites(genome_build),
                   "--fasta", genome_build.reference_fasta, vcf_filename]
    return_code, _stdout, stderr = execute_cmd(extract_cmd)
    if return_code != 0:
        raise ValueError(f"somalier extract failed: {stderr}")

    somalier_files = [os.path.join(extract_dir, f"{sample}.somalier")
                      for sample in (ALLELE_ORDER_SAMPLE, ALLELE_ORDER_UNKNOWN_SAMPLE)]
    relate_cmd = [somalier_bin, "relate", *relate_sites_args, "--unknown", *somalier_files]
    # An unrelated pair is otherwise left out of the report
    env = os.environ | {"SOMALIER_REPORT_ALL_PAIRS": "1"}
    return_code, _stdout, stderr = execute_cmd(relate_cmd, cwd=relate_dir, env=env)
    if return_code != 0:
        raise ValueError(f"somalier relate failed: {stderr}")

    with open(os.path.join(relate_dir, "somalier.pairs.tsv"), encoding="utf-8") as f:
        return int(next(iter(csv.DictReader(f, delimiter="\t")))["ibs0"])


def _allele_order_ibs0(cfg: SomalierConfig, genome_build: GenomeBuild, sites_filename: str,
                       compensate: bool) -> tuple[int, int]:
    """ (number of sites, how many came back as opposite homozygotes) for the probe written and related
        with or without the compensation """
    with tempfile.TemporaryDirectory() as work_dir:
        vcf_filename = os.path.join(work_dir, "allele_order.vcf.bgz")
        num_sites = _write_allele_order_vcf(genome_build, sites_filename, vcf_filename, compensate)
        return_code, _stdout, stderr = execute_cmd(["tabix", vcf_filename])
        if return_code != 0:
            raise ValueError(f"tabix failed: {stderr}")
        return num_sites, _somalier_ibs0(cfg, genome_build, work_dir, vcf_filename, compensate)


def verify_somalier_allele_order() -> dict:
    """ somalier keeps each site's genotypes against its alphabetically sorted alleles rather than the
        record's REF/ALT. v0.3.5's relate --sites only corrected the hom-ref/hom-alt counts in
        somalier.samples.tsv; brentp/somalier edf6c917 (unreleased - we build and host it) also makes
        --unknown impute a no-call as hom-ref the right way round. Consistent relabelling cancels out
        pairwise, so it's only that no-call that shows which somalier is installed: relate a hom-alt
        sample against a no-call one through the installed binary and count opposite homozygotes.

        Uncompensated first: a somalier that gets every site right that way must not be compensated,
        as ancestry compares our extracts with the 1kg ones extracted as called, and the compensation
        turns ours inside out (everyone comes back AMR). Only a somalier that gets it wrong needs
        snpdb.variants_to_vcf to write the AD pair in the site's order, so check that works instead.
        @see snpdb/variants_to_vcf.py:somalier_alleles_flipped """
    cfg = SomalierConfig()
    for genome_build in GenomeBuild.builds_with_annotation():
        sites_filename = cfg.get_sites(genome_build)
        if os.path.exists(sites_filename) and os.path.exists(genome_build.reference_fasta or ""):
            break
    else:
        return {"valid": True,
                "warning": "Somalier allele order unchecked: no build has both a sites file and a fasta"}

    compensated_ibs0 = None
    try:
        try:
            num_sites, uncompensated_ibs0 = _allele_order_ibs0(cfg, genome_build, sites_filename, compensate=False)
        except ValueError as e:
            if RELATE_NO_SITES_OPTION not in str(e):
                raise
            uncompensated_ibs0 = None  # before v0.3.5 - there's no getting it right uncompensated
        if uncompensated_ibs0 is None or uncompensated_ibs0 != num_sites:
            num_sites, compensated_ibs0 = _allele_order_ibs0(cfg, genome_build, sites_filename, compensate=True)
    except Exception as e:
        log_traceback()
        return {"valid": True, "warning": f"Somalier allele order unchecked: {e}"}

    compensating = settings.SOMALIER["compensate_allele_order"]
    return allele_order_result(str(genome_build), num_sites, uncompensated_ibs0, compensated_ibs0, compensating)


def allele_order_result(genome_build_name: str, num_sites: int, uncompensated_ibs0: Optional[int],
                        compensated_ibs0: Optional[int], compensating: bool) -> dict:
    """ Every site is hom-alt against a no-call that --unknown makes hom-ref, so all of them opposite
        is right. uncompensated_ibs0 is None for a somalier with no relate --sites, compensated_ibs0
        None when the uncompensated probe was already right and it wasn't run. """
    needs_compensation = uncompensated_ibs0 != num_sites
    if needs_compensation and compensated_ibs0 != num_sites:
        return {"valid": False,
                "fix": f"{num_sites} hom-alt calls against no-calls came back as {uncompensated_ibs0} "
                       f"(uncompensated) / {compensated_ibs0} (compensated) opposite homozygotes on "
                       f"{genome_build_name}, so somalier is genotyping some third way - check "
                       f"{SOMALIER_ALLELE_ORDER_ISSUE} before trusting relatedness"}
    if needs_compensation == compensating:
        return {"valid": True, "fix": ""}

    if needs_compensation:
        reason = "this somalier reads alleles in alphabetical order, so relatedness is inflated"
    else:
        reason = ("this somalier reads alleles the right way round itself, and compensating turns "
                  "ancestry inside out")
    return {"valid": False,
            "fix": f"{reason} ({genome_build_name}, {SOMALIER_ALLELE_ORDER_ISSUE}). Set "
                   f"SOMALIER[\"compensate_allele_order\"] = {needs_compensation} in this deployment's "
                   "settings and re-run 'somalier_existing_vcfs --clear' to rebuild the extracts"}


def check_somalier() -> dict:
    somalier_data = {
        "somalier_config": {
            "valid": verify_somalier_config(),
            "fix": "Install Somalier, and place in path https://github.com/brentp/somalier/"
        }
    }
    somalier_cfg = SomalierConfig()
    for genome_build in GenomeBuild.builds_with_annotation():
        sites_path = somalier_cfg.get_sites(genome_build)
        somalier_data[f"somalier_sites_file_{genome_build.name}"] = {
            "valid": os.path.exists(sites_path),
            "fix": f"Download somalier sites VCF to {sites_path} - see scripts/install/get_somalier_release.sh",
        }

        try:
            somalier_cfg.get_sites_vcf(genome_build)
            valid = True
        except Exception:
            valid = False
        sites_vcf_name = somalier_cfg.get_sites_vcf_name(genome_build)
        somalier_data[f"somalier_sites_vcf_{genome_build.name}"] = {
            "valid": valid,
            "fix": f"Upload/import '{sites_vcf_name}' as a VCF for {genome_build} (the VCF name must match the sites filename)",
        }
    somalier_data["somalier_allele_order"] = verify_somalier_allele_order()
    return somalier_data
