#!/usr/bin/env python3
"""
We want to do this per-chrom so we can process in parallel

Steps are:
    1. Download exomes.vcf + genome.vcf, removing most INFO fields before writing to disk (to reduce disk space)
    2. Merge exome + genome, summing counts
    3. Run through this script with --af to calculate allele frequency, write TSV (more efficient than VCF)
    4. Cat them all together again
"""

from argparse import ArgumentParser

GRCh38 = "GRCh38"
BUILDS = [GRCh38]

COUNTS = ['AC', 'AN', 'AF']
OTHER_INFOS = ["AC_popmax", "AN_popmax", "AF_popmax", "popmax", "nhomalt", "nhomalt_popmax", "nonpar"]
GNOMAD_SUB_POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]  # Will get AF for each
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'Y']


def get_args():
    parser = ArgumentParser(description="Get, strip and merge gnomAD VCFs for VariantGrid VEP pipeline")
    parser.add_argument("--test", action='store_true', help="Only download 5k of each file.")
    parser.add_argument("--genome-fasta", required=True, help='Fasta (correct for build)')
    parser.add_argument("--chrom_mapping_file", help="Mapping file to convert chroms (if you get 'the sequence 'chr1' was not found)'")
    return parser.parse_args()


def main(args):
    genome_build = GRCh38
    genome_fasta = args.genome_fasta
    if args.test:
        # only download 5k lines of file
        extra_filters = "| bgzip -d | head -5000 | bcftools view -O z"
    else:
        extra_filters = ""  # nothing

    # To remove all INFO tags except "FOO" and "BAR", use "^INFO/FOO,INFO/BAR"
    # @see https://samtools.github.io/bcftools/bcftools.html#annotate """
    info_columns = [f"INFO/{i}" for i in get_columns()]
    keep_columns = ','.join(info_columns)  # AC/AN are special format fields
    bash_header = "#!/bin/bash\nset -e # fail on error\n"

    chrom_scripts = []
    chrom_vcfs = []
    for chrom in CHROMOSOMES:
        prefix = f"gnomad_{genome_build}_chr{chrom}"
        chrom_script = f"{prefix}.sh"
        chrom_scripts.append(chrom_script)
        with open(chrom_script, "w") as cs:
            cs.write(bash_header)

            url = f"https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chr{chrom}.vcf.bgz"
            output_vcf = f"{prefix}.filtered_info.vcf.gz"
            if args.chrom_mapping_file:
                annotate_args = f"--rename-chrs={args.chrom_mapping_file}"
            else:
                annotate_args = ""

            # gnomAD appears to already be decomposed - vt decompose + -s -o +
            cs.write("\necho Download and clean as we go to save disk\n")
            cs.write(f"wget --quiet -O - {url} {extra_filters} | bcftools annotate --exclude 'AC=0' --remove '^{keep_columns}' {annotate_args} | vt normalize - -r {genome_fasta} -o + | vt uniq + -o {output_vcf}\n")

        chrom_vcfs.append(output_vcf)
        if args.test:
            break  # Only do 1 chrom

    # Write merge script
    merge_script_filename = f"gnomad_{genome_build}_merge.sh"


    vcf_header = write_vcf_header()

    with open(merge_script_filename, "w") as ms:
        ms.write(bash_header)
        quoted_files = ' '.join([f"'{f}'" for f in chrom_vcfs])
        gnomad_combined_af_vcf = f"gnomad3_{genome_build}_combined.vcf.bgz"
        ms.write(f"zcat {chrom_vcfs[0]} | head -1000 | grep '^#' | bgzip > vcf_header.bgz")
        ms.write(f"gzcat vcf_header.bgz {quoted_files} | bgzip > {gnomad_combined_af_vcf}\n")
        ms.write(f"tabix {gnomad_combined_af_vcf}\n")

    launch_script_filename = f"gnomad_{genome_build}_launch.sh"
    with open(launch_script_filename, "w") as ms:
        ms.write(bash_header)
        ms.write('SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")\n')
        for cs in chrom_scripts:
            ms.write(f"${{SCRIPT_DIR}}/{cs} > {cs}.log 2> {cs}.stderr.log &\n")

        ms.write("echo Waiting for all chroms to finish...\n")
        ms.write("wait\n")
        ms.write(f"${{SCRIPT_DIR}}/{merge_script_filename}\n")


def write_vcf_header():
    """ Needs to be gzipped so can be concatenated with other gzipped files """
    vcf_header = ""
    return vcf_header


def get_columns():
    columns = COUNTS + OTHER_INFOS
    for g in GNOMAD_SUB_POPS:
        columns.append(f"AF-{g.lower()}")  # gnomAD 3 changed from underscore to dash
    return columns


def write_chrom_mapping_file():
    chrom_mapping_file = "chrom_mapping.txt"
    with open(chrom_mapping_file, "w") as f:
        for c in CHROMOSOMES:
            f.write(f"chr{c}\t{c}\n")
    return chrom_mapping_file


if __name__ == "__main__":
    args = get_args()
    main(args)
