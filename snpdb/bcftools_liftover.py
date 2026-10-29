import logging
import os
from subprocess import CalledProcessError
from typing import Optional

from django.conf import settings

from library.utils import execute_cmd
from snpdb.models import GenomeBuild, VariantCoordinate


def bcftools_lifts_symbolic_as_explicit(variant_coordinate: VariantCoordinate) -> bool:
    """ A del/dup/inv stored symbolic but within LIFTOVER_BCFTOOLS_MAX_LENGTH goes to bcftools as its
        sequence - symbolic storage starts at VARIANT_SYMBOLIC_ALT_SIZE, well below the length bcftools
        lifts explicitly, and its symbolic support is off by default (#1358) """
    if not (variant_coordinate.is_symbolic and variant_coordinate.can_be_made_explicit):
        return False  # explicit already, or <CNV>/<INS> with no sequence to write
    max_length = settings.LIFTOVER_BCFTOOLS_MAX_LENGTH
    return bool(max_length) and abs(variant_coordinate.svlen) < max_length


def bcftools_liftover_variant_coordinate(variant_coordinate: VariantCoordinate, source_genome_build) -> VariantCoordinate:
    """ The form written to the bcftools source VCF. The lifted record is re-imported through the upload
        pipeline, which stores it symbolic again (as_internal_canonical_form) """
    if bcftools_lifts_symbolic_as_explicit(variant_coordinate):
        return variant_coordinate.as_external_explicit(source_genome_build)
    return variant_coordinate


def check_bcftools_coordinate_limits(variant_coordinate: VariantCoordinate) -> Optional[str]:
    """ Returns why the coordinate is outside the settings' symbolic/length limits, or None if bcftools can take it

        bcftools +liftover treats the two forms differently, so their limits differ:
        - Explicit: the REF/ALT sequence is lifted and realigned against the destination fasta, so assembly
          differences inside the variant show up in the result. Cost grows with length (an 11Mb deletion is an
          11Mb REF string), hence LIFTOVER_BCFTOOLS_MAX_LENGTH.
        - Symbolic: lifted like a SNP - only POS is lifted, SVLEN is kept and END is rewritten as POS + the old
          length (freeseek/score#27). The new END is only right when the whole span sits in one gap-free chain
          block, so the error grows with span: a 223Mb <DUP> on chr2 GRCh38->GRCh37 lands END ~870kb away from
          where END itself lifts to. No length check applies here - LIFTOVER_BCFTOOLS_SYMBOLIC is all or nothing.
          Upstream added --lift-end (2025-07, freeseek/score#35) to lift END independently, setting it missing
          when it lands on another contig/strand, but SVLEN is still not updated (#2543) """
    if variant_coordinate.is_symbolic:
        if not settings.LIFTOVER_BCFTOOLS_SYMBOLIC:
            return f"Liftover of symbolic variants disabled via {settings.LIFTOVER_BCFTOOLS_SYMBOLIC=}"
    elif settings.LIFTOVER_BCFTOOLS_MAX_LENGTH:
        if variant_coordinate.max_sequence_length > settings.LIFTOVER_BCFTOOLS_MAX_LENGTH:
            return f"Variant max sequence length > {settings.LIFTOVER_BCFTOOLS_MAX_LENGTH=}"
    return None


def check_bcftools_destination_has_chrom(variant_coordinate: VariantCoordinate,
                                         dest_genome_build: GenomeBuild) -> Optional[str]:
    """ Returns an error if the destination build has no standard contig for the chrom (T2T has no MT: CHM13's
        own mitochondrion is a rotated, non-rCRS sequence). BCFTools would write the record with no header
        contig, failing the entire run """
    if variant_coordinate.chrom in dest_genome_build.chrom_standard_contig_mappings:
        return None
    return f"'{variant_coordinate.chrom}' has no standard contig in {dest_genome_build}"


def check_bcftools_can_liftover(variant_coordinate: VariantCoordinate, source_genome_build: GenomeBuild,
                                dest_genome_build: GenomeBuild) -> Optional[str]:
    """ Returns why bcftools can't liftover the coordinate, or None if it can """
    if error := check_bcftools_coordinate_limits(variant_coordinate):
        return error
    if error := check_bcftools_destination_has_chrom(variant_coordinate, dest_genome_build):
        return error
    # BCFTools requires ref match. Symbolics will pull out reference from build so always match, no point testing
    if not variant_coordinate.is_symbolic:
        calculated_ref = variant_coordinate.calculated_reference(source_genome_build)
        if calculated_ref != variant_coordinate.ref:
            return f"Reference='{variant_coordinate.ref}' not equal to calculated ref from genome {calculated_ref}"
    return None


def bcftools_liftover(source_vcf: str, source_genome_build: GenomeBuild,
                      out_vcf: str, out_genome_build: GenomeBuild) -> tuple[str, int]:
    """ returns reject file and items to process if any failed """

    source_liftover_config = settings.ANNOTATION[source_genome_build.name]["liftover"]
    dest_liftover_config = settings.ANNOTATION[out_genome_build.name]["liftover"]

    source_fasta_filename = source_liftover_config.get("fasta")
    dest_fasta_filename = dest_liftover_config.get("fasta")
    chain_filename = source_liftover_config["chain"].get(out_genome_build.name)
    required_files = {
        "source_fasta_filename": source_fasta_filename,
        "dest_fasta_filename": dest_fasta_filename,
        "chain_filename": chain_filename,
    }

    for description, filename in required_files.items():
        if filename is None:
            raise ValueError(f"BCFTools +liftover couldn't determine what '{description}' to use in liftover from {source_genome_build} -> {out_genome_build}''")
        if not os.path.exists(filename):
            raise ValueError(f"BCFTools +liftover '{description}' ('{filename}') does not exist.")

    reject_vcf = get_reject_vcf_filename(out_vcf)
    bcftools_cmd = "bcftools"  # Installed in path
    cmd = [
        bcftools_cmd,
        "+liftover",
        source_vcf,
        "--",
        "--src-fasta-ref", source_fasta_filename,
        "--fasta-ref", dest_fasta_filename,
        "--chain", chain_filename,
        "--reject", reject_vcf,
        "--write-reject",
    ]

    env = os.environ.copy()
    env["BCFTOOLS_PLUGINS"] = settings.LIFTOVER_BCFTOOLS_PLUGIN_DIR

    logging.info("Executing BCFTools +liftover")
    liftover_cmd_str = " ".join(cmd)
    # bcftools +liftover doesn't sort VCF - @see https://github.com/freeseek/score/issues/8
    sort_cmd = [
        bcftools_cmd,
        "sort",
        "-Oz",
        "-o", out_vcf,
    ]
    sort_cmd_str = " ".join(sort_cmd)
    cmd_str = " | ".join((liftover_cmd_str, sort_cmd_str))
    logging.info(cmd_str)
    # shell=True is required to pipe liftover output into sort. All components are
    # server-side and trusted: the "bcftools" literal plus settings-derived/validated
    # fasta, chain and pipeline VCF paths - no user-supplied input reaches the command.
    return_code, std_out, std_err = execute_cmd([cmd_str], env=env, shell=True)
    logging.info("return_code: %s", return_code)
    if std_out:
        logging.info("stdout: %s", std_out)
    if std_err:
        logging.info("stderr: %s", std_err)

    if return_code:
        raise CalledProcessError(return_code, cmd_str, output=std_out, stderr=std_err)

    items_to_process = 0
    if os.path.exists(reject_vcf):
        items_to_process = count_non_header_lines(reject_vcf)

    return reject_vcf, items_to_process


def get_reject_vcf_filename(out_vcf_filename: str) -> str:
    basename = os.path.basename(out_vcf_filename)
    dirname = os.path.dirname(out_vcf_filename)
    return os.path.join(dirname, "rejected_" + basename)


def count_non_header_lines(filename) -> int:
    count = 0
    with open(filename, encoding="utf-8") as file:
        for line in file:
            if not line.startswith('#'):
                count += 1
    return count
