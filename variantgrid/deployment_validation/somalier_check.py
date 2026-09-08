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
SOMALIER_ALLELE_ORDER_ISSUE = "https://github.com/brentp/somalier/issues/163"


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
    """ One sample, hom-ref at every site, half of them at sites whose ALT sorts first - written the
        way snpdb.variants_to_vcf writes a somalier export. Returns the number of records. """
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

    header_lines = get_vcf_header_from_contigs(genome_build, samples=[ALLELE_ORDER_SAMPLE],
                                               use_accession=False,
                                               formats=[SOMALIER_GT_FORMAT, SOMALIER_AD_FORMAT],
                                               top_lines=[SOMALIER_ALLELE_ORDER_NOTE])
    with open(vcf_filename, "wb") as raw:
        with BGZipWriter(raw) as bgzip_f:
            f = io.TextIOWrapper(bgzip_f, encoding="utf-8", write_through=True)
            writer = VCFWriter(f, header_lines)
            for chrom, position, ref, alt in sites:
                # every call is hom-ref, so all the depth sits on whichever allele we put first
                flipped = somalier_alleles_flipped(ref, alt)
                allele_depths = f"0,{ALLELE_ORDER_DEPTH}" if flipped else f"{ALLELE_ORDER_DEPTH},0"
                writer.write_record(chrom, position, ref, alt, fmt="GT:AD",
                                    sample_calls=[f"0/0:{allele_depths}"])
            f.flush()
            f.detach()

    return len(sites)


def _somalier_genotype_counts(cfg: SomalierConfig, genome_build: GenomeBuild, work_dir: str,
                              vcf_filename: str) -> dict:
    """ extract + relate the one sample, and hand back its row of somalier.samples.tsv """
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

    somalier_file = os.path.join(extract_dir, f"{ALLELE_ORDER_SAMPLE}.somalier")
    return_code, _stdout, stderr = execute_cmd([somalier_bin, "relate", somalier_file], cwd=relate_dir)
    if return_code != 0:
        raise ValueError(f"somalier relate failed: {stderr}")

    with open(os.path.join(relate_dir, "somalier.samples.tsv")) as f:
        return next(iter(csv.DictReader(f, delimiter="\t")))


def verify_somalier_allele_order() -> dict:
    """ somalier reads a record against its own alphabetically sorted site alleles rather than the
        record's REF/ALT, so snpdb.variants_to_vcf writes the AD pair in the site's order to
        compensate. If somalier is ever fixed the compensation becomes the bug and every relatedness
        inverts, and nothing else would say so - a deployment upgrades somalier on its own schedule.
        So genotype a handful of hom-ref calls of each allele order through the installed binary and
        check they come back hom-ref. @see snpdb/variants_to_vcf.py:somalier_alleles_flipped """
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
            counts = _somalier_genotype_counts(cfg, genome_build, work_dir, vcf_filename)
    except Exception as e:
        log_traceback()
        return {"valid": True, "warning": f"Somalier allele order unchecked: {e}"}

    return allele_order_result(str(genome_build), num_sites, int(counts["n_hom_ref"]),
                               int(counts["n_hom_alt"]), settings.SOMALIER["compensate_allele_order"])


def allele_order_result(genome_build_name: str, num_sites: int, hom_ref: int, hom_alt: int,
                        compensating: bool) -> dict:
    """ Every call written was hom-ref, so anything else means the setting and the installed somalier
        disagree about who orders the alleles. Exactly half coming back inverted is that and nothing
        else, and says which way to set it. """
    if hom_ref == num_sites and hom_alt == 0:
        return {"valid": True, "fix": ""}

    fix = (f"{num_sites} hom-ref calls came back as {hom_ref} hom-ref / {hom_alt} hom-alt on "
           f"{genome_build_name} - this somalier and SOMALIER[\"compensate_allele_order\"] "
           f"({compensating}) disagree, see {SOMALIER_ALLELE_ORDER_ISSUE}")
    if hom_ref == ALLELE_ORDER_SITES_PER_ORDER and hom_alt == ALLELE_ORDER_SITES_PER_ORDER:
        fix += (f". Set SOMALIER[\"compensate_allele_order\"] = {not compensating} in this deployment's "
                "settings and re-run 'somalier_existing_vcfs --clear' to rebuild the extracts")
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
