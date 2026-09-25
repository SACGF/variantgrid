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
    somalier_alleles_flipped,
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


def _write_allele_order_vcf(genome_build: GenomeBuild, sites_filename: str, vcf_filename: str) -> int:
    """ One sample hom-alt at every site and one with no call, half of the sites ones whose ALT sorts
        first - written the way snpdb.variants_to_vcf writes a somalier export. Returns the number of
        records. """
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
                flipped = somalier_alleles_flipped(ref, alt)
                allele_depths = f"{ALLELE_ORDER_DEPTH},0" if flipped else f"0,{ALLELE_ORDER_DEPTH}"
                writer.write_record(chrom, position, ref, alt, fmt="GT:AD",
                                    sample_calls=[f"1/1:{allele_depths}", "./.:."])
            f.flush()
            f.detach()

    return len(sites)


def _somalier_ibs0(cfg: SomalierConfig, genome_build: GenomeBuild, work_dir: str, vcf_filename: str) -> int:
    """ extract + relate the pair with --unknown, as the all-samples relate runs, and hand back how many
        sites somalier saw them as opposite homozygotes """
    somalier_bin = cfg.get_annotation("command")
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
    relate_cmd = [somalier_bin, "relate", *cfg.get_relate_sites_args(genome_build), "--unknown", *somalier_files]
    # An unrelated pair is otherwise left out of the report
    env = os.environ | {"SOMALIER_REPORT_ALL_PAIRS": "1"}
    return_code, _stdout, stderr = execute_cmd(relate_cmd, cwd=relate_dir, env=env)
    if return_code != 0:
        raise ValueError(f"somalier relate failed: {stderr}")

    with open(os.path.join(relate_dir, "somalier.pairs.tsv")) as f:
        return int(next(iter(csv.DictReader(f, delimiter="\t")))["ibs0"])


def verify_somalier_allele_order() -> dict:
    """ somalier (up to at least v0.3.5) keeps each site's genotypes against its alphabetically sorted
        alleles rather than the record's REF/ALT, so snpdb.variants_to_vcf writes the AD pair in the
        site's order to compensate. v0.3.5's relate --sites only corrects the hom-ref/hom-alt counts in
        somalier.samples.tsv; pairs are still computed in the sorted order, where --unknown makes a
        missing call hom for whichever allele sorts first. Consistent relabelling cancels out pairwise,
        so it's only that no-call that shows the setting is wrong - and relatedness inflates, with
        nothing else saying so. So relate a hom-alt sample against a no-call one through the installed
        binary, the way the setting says to, and check every site comes back as opposite homozygotes.
        @see snpdb/variants_to_vcf.py:somalier_alleles_flipped """
    cfg = SomalierConfig()
    for genome_build in GenomeBuild.builds_with_annotation():
        sites_filename = cfg.get_sites(genome_build)
        if os.path.exists(sites_filename) and os.path.exists(genome_build.reference_fasta or ""):
            break
    else:
        return {"valid": True,
                "warning": "Somalier allele order unchecked: no build has both a sites file and a fasta"}

    try:
        with tempfile.TemporaryDirectory() as work_dir:
            vcf_filename = os.path.join(work_dir, "allele_order.vcf.bgz")
            num_sites = _write_allele_order_vcf(genome_build, sites_filename, vcf_filename)
            return_code, _stdout, stderr = execute_cmd(["tabix", vcf_filename])
            if return_code != 0:
                raise ValueError(f"tabix failed: {stderr}")
            ibs0 = _somalier_ibs0(cfg, genome_build, work_dir, vcf_filename)
    except Exception as e:
        log_traceback()
        return allele_order_unchecked(str(e))

    compensating = settings.SOMALIER["compensate_allele_order"]
    return allele_order_result(str(genome_build), num_sites, ibs0, compensating)


def allele_order_unchecked(error: str) -> dict:
    """ The probe couldn't run, which is usually a local problem and not something to fail a deploy
        over - except when relate rejected the --sites this deployment's settings say to pass it,
        which is the answer rather than the absence of one. """
    if RELATE_NO_SITES_OPTION in error:
        return {"valid": False,
                "fix": "This somalier predates v0.3.5, whose relate --sites is the only thing that makes "
                       "its hom-ref/hom-alt counts right. Set SOMALIER[\"compensate_allele_order\"] = True "
                       "and re-run 'somalier_existing_vcfs --clear', or install a later somalier. "
                       f"See {SOMALIER_ALLELE_ORDER_ISSUE}"}
    return {"valid": True, "warning": f"Somalier allele order unchecked: {error}"}


def allele_order_result(genome_build_name: str, num_sites: int, ibs0: int, compensating: bool) -> dict:
    """ Every site is hom-alt against a no-call that --unknown makes hom-ref, so anything short of all
        of them opposite means the setting and the installed somalier disagree about who orders the
        alleles. Exactly the half whose ALT sorts first dropping out is that and nothing else, and
        says which way to set it. """
    if ibs0 == num_sites:
        return {"valid": True, "fix": ""}

    fix = (f"{num_sites} hom-alt calls against no-calls came back as {ibs0} opposite homozygotes on "
           f"{genome_build_name} - this somalier and SOMALIER[\"compensate_allele_order\"] "
           f"({compensating}) disagree, see {SOMALIER_ALLELE_ORDER_ISSUE}")
    if ibs0 == ALLELE_ORDER_SITES_PER_ORDER:
        fix += (f". Set SOMALIER[\"compensate_allele_order\"] = {not compensating} in this "
                "deployment's settings and re-run 'somalier_existing_vcfs --clear' to rebuild the extracts")
    else:
        fix += (". Neither setting explains this, so somalier is genotyping some third way - check that "
                "issue before trusting relatedness")
    return {"valid": False, "fix": fix}


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
