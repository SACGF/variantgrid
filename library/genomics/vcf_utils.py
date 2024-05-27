import operator
import os
import re
from collections import defaultdict

import cyvcf2
import vcf

from snpdb.models import Variant, Sequence, GenomeFasta


def cyvcf2_header_types(cyvcf2_reader):
    header_types = defaultdict(dict)
    for h in cyvcf2_reader.header_iter():
        info = h.info()
        h_id = info.get("ID")
        if h_id:  # Not much use w/o this
            header_types[h.type][h_id] = info
    return header_types


def cyvcf2_header_get(cyvcf2_reader, key, default=None):
    try:
        header_dict = cyvcf2_reader.get_header_type(key)
    except:
        header_dict = {}
    return header_dict.get(key, default)


def cyvcf2_get_contig_lengths_dict(cyvcf2_reader):
    try:
        return dict(zip(cyvcf2_reader.seqnames, cyvcf2_reader.seqlens))
    except AttributeError:
        return {}


def write_vcf_from_tuples(vcf_filename, variant_tuples, tuples_have_id_field=False, header_lines: list[str] = None):
    """ variant_tuples can have either 4 or 5 fields (tuples_have_id_field) """

    if header_lines is None:
        header_lines = []

    if tuples_have_id_field:
        vcf_tuples = variant_tuples
    else:
        vcf_tuples = ((chrom, position, ".", ref, alt, svlen) for (chrom, position, ref, alt, svlen) in variant_tuples)

    vcf_tuples = sorted(vcf_tuples, key=operator.itemgetter(0, 1, 3, 4))

    info = [
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    ]
    columns = "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    header = "\n".join(["##fileformat=VCFv4.1", "##source=VariantGrid"] + info + header_lines + ["#" + columns])
    empty_qual_filter_info = ('.', '.', '.')
    with open(vcf_filename, "wt", encoding="utf-8") as f:
        f.write(header + "\n")
        for vcf_record in vcf_tuples:
            (chrom, position, id_col, ref, alt, svlen) = vcf_record
            if Sequence.allele_is_symbolic(alt):
                svtype = alt[1:-1]  # Strip off brackets
                qual_filter_info = (".", ".", f"SVLEN={svlen};SVTYPE={svtype}")
            else:
                qual_filter_info = empty_qual_filter_info
                if alt == Variant.REFERENCE_ALT:
                    alt = "."
            line = "\t".join((chrom, str(position), str(id_col), ref, alt) + qual_filter_info)
            f.write(line + "\n")


def get_variant_caller_and_version_from_vcf(filename) -> tuple[str, str]:
    variant_caller = None
    version = None

    if os.path.exists(filename):
        reader = vcf.Reader(filename=filename)

        if source_list := reader.metadata.get("source"):
            for source in source_list:
                # Match source = "freeBayes v1.3.5" or "VarDict_v1.8.2"
                if m := re.match(r"(.*?)[ _]v([\d\\.]+)", source):
                    variant_caller, version = m.groups()
                    break

        if gatk_commandline := reader.metadata.get("GATKCommandLine"):
            variant_caller = "GATK"
            for commandline in gatk_commandline:
                if caller_id := commandline.get("ID"):
                    if caller_id == 'HaplotypeCaller':
                        caller_id = "GATK"  # Just stay with GATK
                    variant_caller = caller_id

                if version := commandline.get("Version"):
                    version = version.replace('"', "")  # Strip quotes

    return variant_caller, version


def vcf_allele_is_symbolic(allele: str) -> bool:
    return allele.startswith("<") and allele.endswith(">")


def vcf_get_ref_alt_svlen(variant: cyvcf2.Variant):
    ref = variant.REF.strip().upper()
    if variant.ALT:
        alt = variant.ALT[0].strip().upper()
    else:
        alt = Variant.REFERENCE_ALT

    if Sequence.allele_is_symbolic(ref) or Sequence.allele_is_symbolic(alt):
        # Need to provide END or SVLEN
        if svlen_info := variant.INFO.get('SVLEN'):
            svlen = int(svlen_info)
        elif end_info := variant.INFO.get('END'):
            svlen = int(end_info) - variant.POS
        else:
            raise ValueError(f"SVLEN or END info field MUST be provided for symbolic (ie '<x>') {ref=},{alt=}")
    else:
        svlen = None
    return ref, alt, svlen


def get_vcf_header_contig_lines(contigs: list[tuple]) -> list[str]:
    header_lines = []
    for contig, length, assembly in contigs:
        line = f"##contig=<ID={contig},length={length},assembly={assembly}>"
        header_lines.append(line)
    return header_lines


def get_contigs_header_lines(genome_build, standard_only=True, contig_allow_list: set = None) -> list[str]:
    if standard_only:
        contig_qs = genome_build.standard_contigs
    else:
        contig_qs = genome_build.contigs

    contigs = []
    for contig in contig_qs:
        if contig_allow_list is not None:
            if contig.refseq_accession not in contig_allow_list:
                continue
        contigs.append((contig.refseq_accession, contig.length, genome_build.name))
    return get_vcf_header_contig_lines(contigs)


def write_cleaned_vcf_header(genome_build, source_vcf_filename: str, output_filename: str, new_info_lines: list[str] = None):
    contig_regex = re.compile(r"^##contig=<ID=(.+),length=(\d+)")

    header_lines = []
    with open(source_vcf_filename) as in_f:
        for line in in_f:
            if not line.startswith("#"):
                break  # End of header
            header_lines.append(line.strip())

    # These are used to validate contigs in header
    genome_fasta = GenomeFasta.get_for_genome_build(genome_build)
    chrom_to_contig_id = genome_build.get_chrom_contig_id_mappings()
    contig_lengths = dict(genome_build.contigs.values_list("pk", "length"))
    contig_to_fasta_names = genome_fasta.get_contig_id_to_name_mappings()

    with open(output_filename, "w") as f:
        found_column_names_line = False
        for line in header_lines:
            if line.startswith("#CHROM"):
                found_column_names_line = True
                # This is where we dump the new stuff
                if new_info_lines:
                    for new_info_line in new_info_lines:
                        f.write(new_info_line + "\n")
                for contig_line in get_contigs_header_lines(genome_build):
                    f.write(contig_line + "\n")

            elif m := contig_regex.match(line):
                # Strip existing contig lines from header - though check they match so we don't get build swaps
                contig_name, provided_contig_length = m.groups()
                if contig_id := chrom_to_contig_id.get(contig_name):
                    if fasta_chrom := contig_to_fasta_names.get(contig_id):
                        provided_contig_length = int(provided_contig_length)
                        ref_contig_length = contig_lengths[contig_id]
                        if provided_contig_length != ref_contig_length:
                            msg = f"VCF header contig '{contig_name}' (length={provided_contig_length}) has " + \
                                f"different length than ref contig {fasta_chrom} (length={ref_contig_length})"
                            raise ValueError(msg)

            f.write(line + "\n")

        if not found_column_names_line:
            raise ValueError("VCF header was missing line starting with '#CHROM'")
