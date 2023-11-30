#!/usr/bin/env python3
"""
We want to do this per-chrom so we can process in parallel

Steps are:
    1. Download exomes.vcf + genome.vcf, removing most INFO fields before writing to disk (to reduce disk space)
    2. Merge exome + genome, summing counts
    3. Run through this script with --af to calculate allele frequency, write TSV (more efficient than VCF)
    4. Cat them all together again
"""

import gzip
import os
from argparse import ArgumentParser
from datetime import datetime
from typing import Tuple, List

GNOMAD_V_2_1 = "2.1.1"
GNOMAD_V_3_1_2 = "3.1.2"
GNOMAD_V_4_0 = "4.0"

GNOMAD_VERSIONS = {
    GNOMAD_V_2_1,
    GNOMAD_V_3_1_2,
    GNOMAD_V_4_0,
}

FILENAMES = {
    GNOMAD_V_2_1: "gnomad.%(capture_type)s.r2.1.1.sites.%(chrom)s.vcf.bgz",
    GNOMAD_V_3_1_2: "gnomad.%(capture_type)s.v3.1.2.sites.chr{chrom}.vcf.bgz",
    GNOMAD_V_4_0: "gnomad.%(capture_type)s.v4.0.sites.chr%(chrom)s.vcf.bgz",
}


GENOME_BUILDS = {"GRCh37", "GRCh38"}

def get_infos_for_version(gnomad_version) -> Tuple[List[str], List[str], List[str], List[str]]:
    # We deliberately leave out AF and "grpmax" stuff as we recalculate that later in 'calculate_allele_frequency'
    info_fields = ['AC', 'AN', "nhomalt", "nonpar"]
    popmax_fields = ["AF_popmax", "AC_popmax", "AN_popmax", "popmax", "nhomalt_popmax"]
    grpmax_fields = ["AF_grpmax", "AC_grpmax", "AN_grpmax", "grpmax", "nhomalt_grpmax"]
    sub_pops = ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]  # Will get AF for each
    chr_x_male = ["AC_male", "AN_male", "AF_male"]
    chr_x_xy = ["AC_XY", "AN_XY", "AF_XY"]

    if gnomad_version == GNOMAD_V_4_0:
        popmax_fields = grpmax_fields
        chr_x_male = chr_x_xy
        info_fields.extend(["faf95", "faf99", "fafmax_faf95_max", "fafmax_faf99_max"])
        # Others are now called remaining
        sub_pops.remove("oth")
        sub_pops.append("remaining")  #
        sub_pops.append("mid")  # Middle easterners added in v4

        info_fields.remove("nonpar")
        info_fields.append("non_par")

    return info_fields, chr_x_male, popmax_fields, sub_pops

# popmax/grpmax is calculated using non-bottlenecked genetic ancestry groups
BOTTLENECKED_SUB_POPS = {"asj", "fin", "mid", "oth", "remaining"}


def get_args():
    available_builds = ", ".join(GENOME_BUILDS)
    available_versions = ", ".join(GNOMAD_VERSIONS)

    parser = ArgumentParser(description="Merge exome+genome VCFs for VariantGrid VEP pipeline")
    parser.add_argument("--test", action='store_true', help="Only do chrY (quick test)")
    parser.add_argument("--chrom-mapping-file", help='bcftools chromosome conversion')
    parser.add_argument("--genome-build", help=f'GenomeBuild (one of {available_builds})')
    parser.add_argument("--version", help=f'gnomAD version (one of {available_versions})')
    parser.add_argument("--path", help='Optional Colon separated paths for tabix/bgzip/vt/bcftools')
    parser.add_argument("--gnomad-input-vcf")
    parser.add_argument("--af-output-vcf")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--scripts', action='store_true', help="Generate scripts")
    group.add_argument('--af', action='store_true', help="Calculate allele frequency from VCF")

    args = parser.parse_args()
    if args.scripts:
        if args.genome_build is None:
            parser.error("--genome-build required for --scripts")
        if args.genome_build not in GENOME_BUILDS:
            parser.error(f"--genome-build must be one of {', '.join(GENOME_BUILDS)}")
    else:
        if args.gnomad_input_vcf is None:
            parser.error("--gnomad-input-vcf required for --af")
        if args.af_output_vcf is None:
            parser.error("--af-output-vcf required for --af")

    if args.version not in GNOMAD_VERSIONS:
        parser.error(f"Version must be one of: {available_versions}")

    return args


def main(args):
    if args.scripts:
        write_scripts(args)
    else:
        calculate_allele_frequency(args.version, args.gnomad_input_vcf, args.af_output_vcf)


def write_scripts(args):
    if not args.chrom_mapping_file:
        raise ValueError("--chrom-mapping-file is required for write scripts step")

    if args.test:
        chromosomes = ["Y"]  # Just do Y
    else:
        chromosomes = list(map(str, range(1, 23))) + ['X', 'Y']

    info_fields, chr_x_male, popmax_fields, sub_pops = get_infos_for_version(args.version)

    columns = get_columns(info_fields, sub_pops)
    bash_header = "#!/bin/bash\nset -e # fail on error\n"

    if args.path:
        bash_header += "PATH=${PATH}:" + args.path + "\n"

    filename_template = FILENAMES[args.version]

    chrom_scripts = []
    af_vcfs = []
    for chrom in chromosomes:
        prefix = f"gnomad{args.version}_{args.genome_build}_chr{chrom}"
        chrom_script = f"{prefix}.sh"
        chrom_scripts.append(chrom_script)
        with (open(chrom_script, "w") as cs):
            cs.write(bash_header)

            output_vcfs = []
            for capture_type in ["exomes", "genomes"]:
                # To remove all INFO tags except "FOO" and "BAR", use "^INFO/FOO,INFO/BAR"
                # @see https://samtools.github.io/bcftools/bcftools.html#annotate """
                my_columns = columns.copy()
                if chrom == "X":
                    my_columns.extend(chr_x_male)

                info_columns = [f"INFO/{i}" for i in my_columns]
                keep_columns = ','.join(info_columns)  # AC/AN are special format fields
                output_vcf = f"{prefix}_{capture_type}.filtered_info.vcf.gz"
                annotate_args = f"--rename-chrs={args.chrom_mapping_file}"

                gnomad_vcf_filename = filename_template % {
                    "capture_type": capture_type,
                    "chrom": chrom,
                }

                # bcftools now works with AC/AN etc - see https://github.com/samtools/bcftools/issues/1394
                # but make sure you are using v18
                # bcftools merge doesn't work with type='A'
                modify_fields2 = "sed -e 's/,Number=A,/,Number=1,/'"
                # gnomAD appears to already be decomposed - vt decompose + -s -o +
                # We no longer remove AC=0 as we want to keep AN (total counts) for pops for later AF calculations
                cs.write(f"bcftools annotate --remove '^{keep_columns}' {annotate_args} {gnomad_vcf_filename} | {modify_fields2} | vt uniq + -o {output_vcf}\n")
                output_vcfs.append(output_vcf)

            combined_vcf = f"{prefix}.combined.vcf.gz"
            if len(output_vcfs) == 1:  # Just 1, rename it
                output_vcf = output_vcfs[0]
                cs.write(f"mv {output_vcf} {combined_vcf}\n")
            else:
                for ov in output_vcfs:
                    cs.write(f"tabix {ov}\n")

                # Merge exomes/genome VCFs
                # if we leave out rule, will take from 1st file which is ok for PAR as will be the same
                skip_columns = {"nonpar", "non_par"}
                # Default rule = "sum" if not below (or skipped)
                rule_ops = {
                    # Will take higher of whatever is there in genomes/exomes
                    "faf95": "max",
                    "faf99": "max",
                    "fafmax_faf95_max": "max",
                    "fafmax_faf99_max": "max",
                }
                info_rules = []
                for c in my_columns:
                    if c not in skip_columns:
                        op = rule_ops.get(c, "sum")
                        info_rules.append(f"{c}:{op}")
                info_rules_arg = ','.join(info_rules)
                cs.write("\n\necho Merging VCFs - will keep flags from genomes.\n")
                # Throw away anything that has AC=0
                cs.write(f"bcftools merge --merge none --info-rules '{info_rules_arg}' '{output_vcfs[0]}' '{output_vcfs[1]}' | bcftools view -i 'AC>0' -O z -o {combined_vcf}\n")

            # Now process them with this script
            cs.write("\n\necho Calculate Allele Frequency\n")
            # Lines for Phoenix HPC
            # cs.write("\nmodule load Python/3.9.6-GCCcore-11.2.0\n")
            # cs.write("source /home/a1059391/venv/dave_venv/bin/activate\n")
            script_filename = os.path.realpath(__file__)
            allele_frequency_vcf = f"{prefix}.af.vcf.gz"
            cs.write(f"{script_filename} --af --gnomad-input-vcf={combined_vcf} --af-output-vcf={allele_frequency_vcf} --version={args.version}\n")
            af_vcfs.append(allele_frequency_vcf)

    # Write merge script
    merge_script_filename = f"gnomad{args.version}_merge.sh"
    vcf_header = write_vcf_header(args.version, info_fields, popmax_fields, sub_pops)

    with open(merge_script_filename, "w") as ms:
        ms.write(bash_header)
        quoted_files = ' '.join([f"'{f}'" for f in af_vcfs])
        gnomad_combined_af_vcf = f"gnomad{args.version}_{args.genome_build}_combined_af.vcf.bgz"
        # We produce gzipped files, but want bgzipped, so need to cat then bgzip
        ms.write(f"zcat {vcf_header} {quoted_files} | bgzip > {gnomad_combined_af_vcf}\n")
        ms.write(f"tabix {gnomad_combined_af_vcf}\n")

    launch_script_filename = f"gnomad{args.version}_launch.sh"
    with open(launch_script_filename, "w") as ms:
        ms.write(bash_header)
        ms.write('SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")\n')
        for cs in chrom_scripts:
            ms.write(f"${{SCRIPT_DIR}}/{cs} > {cs}.log 2> {cs}.stderr.log &\n")

        ms.write("echo Waiting for all chroms to finish...\n")
        ms.write("wait\n")
        ms.write(f"${{SCRIPT_DIR}}/{merge_script_filename}\n")


def get_columns(info_fields, sub_pops):
    columns = info_fields.copy()
    for g in sub_pops:
        for f in ["AC", "AN"]:
            columns.append(f"{f}_{g.lower()}")
    return columns


def get_af_info(sub_pops):
    af_info = [
        ("AF", None, "AC", "AN"),
    ]
    for g in sub_pops:
        af_info.append((f'AF_{g}', g, f'AC_{g}', f'AN_{g}'))
    return af_info


def write_vcf_header(version, info_fields, popmax_fields, sub_pops):
    """ Needs to be gzipped so can be concatenated with other gzipped files """

    all_info = set(info_fields + popmax_fields + ["gnomad_filtered"])
    field_headers = {
        'AC': '##INFO=<ID=AC,Number=1,Type=Integer,Description="Alternate allele count (exomes + genomes)">',
        'AN': '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles  (exomes + genomes)">',
        'AC_XY': '##INFO=<ID=AC_XY,Number=1,Type=Integer,Description="Alternate allele count for XY samples">',
        'AF_XY': '##INFO=<ID=AF_XY,Number=1,Type=Float,Description="Alternate allele frequency in XY samples">',
        'AN_XY': '##INFO=<ID=AN_XY,Number=1,Type=Integer,Description="Total number of alleles in XY samples">',
        'AC_male': '##INFO=<ID=AC_male,Number=1,Type=Integer,Description="Alternate allele count for male samples">',
        'AN_male': '##INFO=<ID=AN_male,Number=1,Type=Integer,Description="Total number of alleles in male samples">',
        'AF_male': '##INFO=<ID=AF_male,Number=1,Type=Float,Description="Alternate allele frequency in male samples">',
        'faf95': '##INFO=<ID=faf95,Number=1,Type=Float,Description="Filtering allele frequency (using Poisson 95% CI) (max of exomes/genomes)">',
        'faf99': '##INFO=<ID=faf99,Number=1,Type=Float,Description="Filtering allele frequency (using Poisson 99% CI) (max of exomes/genomes)">',
        'fafmax_faf95_max': '##INFO=<ID=fafmax_faf95_max,Number=1,Type=Float,Description="Maximum filtering allele frequency (using Poisson 95% CI) across genetic_ancestry groups (max of exomes/genomes)">',
        'fafmax_faf99_max': '##INFO=<ID=fafmax_faf99_max,Number=1,Type=Float,Description="Maximum filtering allele frequency (using Poisson 99% CI) across genetic_ancestry groups (max of exomes/genomes)">',
        'AF_popmax': '##INFO=<ID=AF_popmax,Number=1,Type=Float,Description="Allele Frequency for highest population">',
        'AC_popmax': '##INFO=<ID=AC_popmax,Number=1,Type=Integer,Description="Allele Count for highest population">',
        'AN_popmax': '##INFO=<ID=AN_popmax,Number=1,Type=Integer,Description="Allele Number for highest population">',
        'popmax': '##INFO=<ID=popmax,Number=1,Type=String,Description="Ancestral group with highest allele frequency (stored as AF_popmax)">',
        'AF_grpmax': '##INFO=<ID=AF_grpmax,Number=1,Type=Float,Description="Allele Frequency for highest population">"',
        'AC_grpmax': '##INFO=<ID=AC_grpmax,Number=1,Type=Integer,Description="Allele Count for highest population">',
        'AN_grpmax': '##INFO=<ID=AN_grpmax,Number=1,Type=Integer,Description="Allele Number for highest population">',
        'grpmax': '##INFO=<ID=grpmax,Number=1,Type=String,Description="Ancestral group with highest allele frequency (stored as AF_grpmax)">',
        'nhomalt': '##INFO=<ID=nhomalt,Number=1,Type=Integer,Description="Total number of homozygotest (exomes + genomes)">',
        'gnomad_filtered': '##INFO=<ID=gnomad_filtered,Number=1,Type=Integer,Description="Exomes or genomes had a filter entry (potential QC issues)">',
        'nonpar': '##INFO=<ID=nonpar,Number=0,Type=Flag,Description="Variant (on sex chromosome) falls outside a pseudoautosomal region">',
        'non_par': '##INFO=<ID=nonpar,Number=0,Type=Flag,Description="Variant (on sex chromosome) falls outside a pseudoautosomal region">',
    }

    info_headers = ""
    for field in all_info:
        if header := field_headers.get(field):
            info_headers += header + "\n"

    now = datetime.now()
    file_date = "%d%02d%02d" % (now.year, now.month, now.day)
    source = __file__
    meta = """##fileformat=VCFv4.2
##fileDate=%(file_date)s
##source=%(source)s
%(info_headers)s""" % {"file_date": file_date, "source": source, "info_headers": info_headers}

    af_info = get_af_info(sub_pops)
    for info_id, pop_name, ac_name, an_name in af_info:
        if pop_name:
            af_desc = f"for {pop_name}"
        else:
            af_desc = ""
        af_desc += f" made from (exomes_{ac_name} + genomes_{ac_name}) / (exomes_{an_name} + genomes_{an_name})"
        meta += f'##INFO=<ID={info_id},Number=1,Type=Float,Description="Allele Frequency {af_desc}">\n'

    vcf_header = f"gnomad_{version}_vcf_header.txt.gz"
    with gzip.open(vcf_header, "wt") as f:
        f.write(meta)
        header_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        header = "#" + '\t'.join(header_cols)
        f.write(header + "\n")
    return vcf_header


def calculate_allele_frequency(version, gnomad_input_vcf, af_output_vcf):
    """ We have to re-calculate POPMAX as we can't merge it """

    from cyvcf2 import VCF  # Import here, so that rest of script can run on HPC easier

    info_fields, chr_x_male, popmax_fields, sub_pops = get_infos_for_version(version)
    af_info = get_af_info(sub_pops)

    with gzip.open(af_output_vcf, "wt") as f:
        for variant in VCF(gnomad_input_vcf):
            chrom = variant.CHROM
            pos = str(variant.POS)
            variant_id = variant.ID or '.'
            ref = variant.REF
            alt = variant.ALT[0]  # no multi-alts

            af_popmax = 0
            ac_popmax = 0
            an_popmax = 0
            popmax = '.'
            infos = {}
            for af_name, pop_name, ac_name, an_name in af_info:
                ac = variant.INFO.get(ac_name, 0)
                an = variant.INFO.get(an_name)
                # print(f"{pop_name=},{ac_name=},{an_name=} {ac=}/{an=}")
                if an:
                    af = ac / an
                    if pop_name and (pop_name not in BOTTLENECKED_SUB_POPS) and af > af_popmax:
                        af_popmax = af
                        ac_popmax = ac
                        an_popmax = an
                        popmax = pop_name
                    af = f'{af:.6f}'
                else:
                    af = '.'

                infos[af_name] = af
                infos[ac_name] = ac
                infos[an_name] = an

            for o in info_fields + chr_x_male:
                infos[o] = str(variant.INFO.get(o, '.'))
            gnomad_filtered = '0' if variant.FILTER is None else '1'
            infos["gnomad_filtered"] = gnomad_filtered

            for p in popmax_fields: # can be popmax or grpmax
                if p.startswith("AF_"):
                    infos[p] = str(af_popmax)
                elif p.startswith("AC_"):
                    infos[p] = str(ac_popmax)
                elif p.startswith(("AN_")):
                    infos[p] = str(an_popmax)
                elif p in {"popmax", "grpmax"}:
                    infos[p] = popmax
            info_str = ";".join([f"{k}={v}" for k, v in infos.items()])
            columns = [chrom, pos, variant_id, ref, alt, '.', '.', info_str]
            f.write("\t".join(columns) + "\n")


if __name__ == "__main__":
    args = get_args()
    main(args)
