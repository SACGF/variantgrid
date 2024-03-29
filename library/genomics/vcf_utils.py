import operator
import os
import re
from collections import defaultdict

import cyvcf2
import vcf

from snpdb.models import Variant, Sequence


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


def write_vcf_from_tuples(vcf_filename, variant_tuples, tuples_have_id_field=False):
    """ variant_tuples can have either 4 or 5 fields (tuples_have_id_field) """

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
    header = "\n".join(["##fileformat=VCFv4.1", "##source=VariantGrid"] + info + ["#" + columns])
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
