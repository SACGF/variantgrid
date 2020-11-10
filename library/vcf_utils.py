import operator
from collections import defaultdict


class VCFConstant:
    FREEBAYES = "freeBayes"
    CLCAD2 = "CLCAD2"  # CLC Genomics Workbench - variant track counts for ref,alt (1 for each alt)
    DEFAULT_ALLELE_FIELD = 'AD'
    DEFAULT_READ_DEPTH_FIELD = 'DP'
    DEFAULT_GENOTYPE_QUALITY_FIELD = 'GQ'
    DEFAULT_PHRED_LIKILIHOOD_FIELD = 'PL'
    GENOTYPE_LIKELIHOOD = "GL"
    ALT_DEPTH_FIELD = "AO"  # FreeBayes - Alternate allele observation count
    REF_DEPTH_FIELD = "RO"  # FreeBayes - Reference allele observation count


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
        vcf_tuples = ((chrom, position, ".", ref, alt) for (chrom, position, ref, alt) in variant_tuples)

    vcf_tuples = sorted(vcf_tuples, key=operator.itemgetter(0, 1, 3, 4))
    columns = "\t".join(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
    header = "\n".join(["##fileformat=VCFv4.1", "##source=VariantGrid", "#" + columns])
    empty_columns = "\t." * 3  # QUAL/FILTER/INFO always empty
    with open(vcf_filename, "wt") as f:
        f.write(header + "\n")
        for vcf_record in vcf_tuples:
            line = "\t".join([str(s) for s in vcf_record]) + empty_columns
            f.write(line + "\n")
