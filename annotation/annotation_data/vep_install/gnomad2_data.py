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
from cyvcf2 import VCF
from datetime import datetime
import gzip
import os

GRCh37 = "GRCh37"
GRCh38 = "GRCh38"
BUILDS = [GRCh37, GRCh38]

COUNTS = ['AC', 'AN']
OTHER_INFOS = ["nhomalt"]
GNOMAD_SUB_POPS = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]  # Will get AN/AC for these
CHROMOSOMES = list(map(str, range(1, 23))) + ['X', 'Y']


def get_args():
    parser = ArgumentParser(description="Merge exome+genome VCFs for VariantGrid VEP pipeline")
    parser.add_argument("--test", action='store_true', help="Only download 5k of each file.")
    parser.add_argument("--genome-build", help='GRCh37 or GRCh38')
    parser.add_argument("--genome-fasta", help='Fasta (correct for build)')
    parser.add_argument("--genome-fasta-has-chr", type=bool, help='GnomAD has "chr1", set to false if ref uses "1"')
    parser.add_argument("--version", help='gnomAD version (default: 2.1.1)', default='2.1.1')
    parser.add_argument("--gnomad-input-vcf")
    parser.add_argument("--af-output-vcf")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--scripts', action='store_true', help="Generate scripts")
    group.add_argument('--af', action='store_true', help="Calculate allele frequency from VCF")

    args = parser.parse_args()
    if args.scripts:
        if args.genome_build is None or args.genome_build not in BUILDS:
            parser.error(f"--genome-build must be one of {','.join(BUILDS)}")
        if args.genome_fasta is None:
            parser.error("--genome-fasta required for --scripts")
    else:
        if args.gnomad_input_vcf is None:
            parser.error("--gnomad-input-vcf required for --af")
        if args.af_output_vcf is None:
            parser.error("--af-output-vcf required for --af")

    return args


def main(args):
    if args.scripts:
        write_scripts(args)
    else:
        calculate_allele_frequency(args.gnomad_input_vcf, args.af_output_vcf)


def write_scripts(args):
    genome_build = args.genome_build
    version = args.version
    genome_fasta = args.genome_fasta
    if args.test:
        # only download 5k lines of file
        extra_filters = "| bgzip -d | head -5000 | bcftools view -O z"
    else:
        extra_filters = ""  # nothing

    if not args.genome_fasta_has_chr:
        chrom_mapping_file = None
    else:
        chrom_mapping_file = write_chrom_mapping_file()

    columns = get_columns()
    bash_header = "#!/bin/bash\nset -e # fail on error\n"

    chrom_scripts = []
    af_vcfs = []
    for chrom in CHROMOSOMES:
        prefix = f"gnomad_{genome_build}_chr{chrom}"
        chrom_script = f"{prefix}.sh"
        chrom_scripts.append(chrom_script)
        with open(chrom_script, "w") as cs:
            cs.write(bash_header)

            output_vcfs = []
            for vcf_type in ["exomes", "genomes"]:
                if genome_build == GRCh37:
                    url = f"https://storage.googleapis.com/gnomad-public/release/{version}/vcf/{vcf_type}/gnomad.{vcf_type}.r{version}.sites.{chrom}.vcf.bgz"
                else:
                    url = f"https://storage.googleapis.com/gnomad-public/release/{version}/liftover_grch38/vcf/{vcf_type}/gnomad.{vcf_type}.r{version}.sites.{chrom}.liftover_grch38.vcf.bgz"

                # To remove all INFO tags except "FOO" and "BAR", use "^INFO/FOO,INFO/BAR"
                # @see https://samtools.github.io/bcftools/bcftools.html#annotate """
                my_columns = columns.copy()
                if vcf_type == "genomes":  # No SAS in genomes
                    if chrom == 'Y':
                        continue  # No Y in genomes
                    my_columns.remove("AC_sas")
                    my_columns.remove("AN_sas")

                info_columns = [f"INFO/{i}" for i in my_columns]
                keep_columns = ','.join(info_columns)  # AC/AN are special format fields
                output_vcf = f"{prefix}_{vcf_type}.filtered_info.vcf.gz"
                if chrom_mapping_file:
                    annotate_args = f"--rename-chrs={chrom_mapping_file}"
                else:
                    annotate_args = ""

                # bcftools merge doesn't work with type='A' or special AC/AN INFO fields w/o a FORMAT (which gnomAD doesn't have)
                modify_fields = "sed -e 's/,Number=A,/,Number=1,/' -e 's/ID=AC,/ID=AC_count,/' -e 's/ID=AN,/ID=AN_count,/' -e 's/AC=/AC_count=/' -e 's/AN=/AN_count=/'"
                # gnomAD appears to already be decomposed - vt decompose + -s -o +
                cs.write("\necho Download and clean as we go to save disk\n")
                cs.write(f"wget --quiet -O - {url} {extra_filters} | bcftools annotate --exclude 'AC=0' --remove '^{keep_columns}' {annotate_args} | {modify_fields} | vt normalize - -r {genome_fasta} -o + | vt uniq + -o {output_vcf}\n")
                output_vcfs.append(output_vcf)

            combined_vcf = f"{prefix}.combined.vcf.gz"
            if len(output_vcfs) == 1:  # Just 1, rename it
                output_vcf = output_vcfs[0]
                cs.write(f"mv {output_vcf} {combined_vcf}\n")
            else:
                for ov in output_vcfs:
                    cs.write(f"tabix {ov}\n")

                # Merge - adding them together...
                renamed_columns = [f"{c}_count" if c in ['AC', 'AN'] else c for c in columns]
                info_rules = [f"{c}:sum" for c in renamed_columns]
                info_rules_arg = ','.join(info_rules)
                cs.write("\n\necho Merging VCFs - will keep flags from genomes.\n")
                cs.write(f"bcftools merge --merge none --info-rules '{info_rules_arg}' '{output_vcfs[0]}' '{output_vcfs[1]}' -O z -o {combined_vcf}\n")

            # Now process them with this script
            cs.write("\n\necho Calculate Allele Frequency\n")
            script_filename = os.path.realpath(__file__)
            allele_frequency_vcf = f"{prefix}.af.vcf.gz"
            cs.write(f"{script_filename} --af --gnomad-input-vcf={combined_vcf} --af-output-vcf={allele_frequency_vcf}\n")
            af_vcfs.append(allele_frequency_vcf)

        if args.test:
            break  # Only do 1 chrom

    # Write merge script
    merge_script_filename = f"gnomad_{genome_build}_merge.sh"
    vcf_header = write_vcf_header()

    with open(merge_script_filename, "w") as ms:
        ms.write(bash_header)
        quoted_files = ' '.join([f"'{f}'" for f in af_vcfs])
        gnomad_combined_af_vcf = f"gnomad_{genome_build}_combined_af.vcf.bgz"
        ms.write(f"gzcat {vcf_header} {quoted_files} | bgzip > {gnomad_combined_af_vcf}\n")
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


def get_columns():
    columns = COUNTS + OTHER_INFOS
    for g in GNOMAD_SUB_POPS:
        for f in ["AC", "AN"]:
            columns.append(f"{f}_{g.lower()}")
    return columns


def get_af_info():
    af_info = [
        ("AF", None, "AC_count", "AN_count"),
    ]
    for g in GNOMAD_SUB_POPS:
        af_info.append((f'AF_{g}', g, f'AC_{g}', f'AN_{g}'))
    return af_info


def write_vcf_header():
    """ Needs to be gzipped so can be concatenated with other gzipped files """

    now = datetime.now()
    file_date = "%d%02d%02d" % (now.year, now.month, now.day)
    source = __file__
    meta = """##fileformat=VCFv4.2
##fileDate=%(file_date)s
##source=%(source)s
##INFO=<ID=AF_popmax,Number=1,Type=Float,Description="Allele Frequency for highest population">
##INFO=<ID=popmax,Number=1,Type=String,Description="Population with highest allele frequency (stored as AF_popmax)">
##INFO=<ID=nhomalt,Number=1,Type=Integer,Description="Total number of homozygotest (exomes + genomes)">
##INFO=<ID=gnomad_filtered,Number=1,Type=Integer,Description="Exomes or genomes had a filter entry (potential QC issues)">
""" % {"file_date": file_date, "source": source}

    af_info = get_af_info()
    for info_id, pop_name, ac_name, an_name in af_info:
        if pop_name:
            af_desc = f"for {pop_name}"
        else:
            af_desc = ""
        af_desc += f" made from (exomes_{ac_name} + genomes_{ac_name}) / (exomes_{an_name} + genomes_{an_name})"
        meta += f'##INFO=<ID={info_id},Number=1,Type=Float,Description="Allele Frequency {af_desc}">\n'

    vcf_header = "vcf_header.txt.gz"
    with gzip.open(vcf_header, "wt") as f:
        f.write(meta)
        header_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        header = "#" + '\t'.join(header_cols)
        f.write(header + "\n")
    return vcf_header


def write_chrom_mapping_file():
    chrom_mapping_file = "chrom_mapping.txt"
    with open(chrom_mapping_file, "w") as f:
        for c in CHROMOSOMES:
            f.write(f"chr{c}\t{c}\n")
    return chrom_mapping_file


def calculate_allele_frequency(gnomad_input_vcf, af_output_vcf):
    # We have to re-calculate POPMAX as we can't merge it
    af_info = get_af_info()
    info_names = [ai[0] for ai in af_info] + OTHER_INFOS + ["AF_popmax", "popmax", "gnomad_filtered"]

    with gzip.open(af_output_vcf, "wt") as f:
        for variant in VCF(gnomad_input_vcf):
            chrom = variant.CHROM
            pos = str(variant.POS)
            variant_id = variant.ID or '.'
            ref = variant.REF
            alt = variant.ALT[0]  # no multi-alts

            af_popmax = 0
            popmax = '.'
            infos = []
            for _, pop_name, ac_name, an_name in af_info:
                ac = variant.INFO.get(ac_name, 0)
                an = variant.INFO.get(an_name)
                #print(f"{ac_name}/{an_name} {ac}/{an}")
                if an:
                    af = ac / an
                    if pop_name and af > af_popmax:  # Only use subpops
                        af_popmax = af
                        popmax = pop_name
                    af = '{:.6f}'.format(af)
                else:
                    af = '.'
                infos.append(af)

            for o in OTHER_INFOS:
                infos.append(str(variant.INFO.get(o, '.')))
            gnomad_filtered = '0' if variant.FILTER is None else '1'
            infos.extend([str(af_popmax), popmax, gnomad_filtered])
            info_str = ";".join([i + "=" + v for i, v in zip(info_names, infos)])
            columns = [chrom, pos, variant_id, ref, alt, '.', '.', info_str]
            f.write("\t".join(columns) + "\n")


if __name__ == "__main__":
    args = get_args()
    main(args)
