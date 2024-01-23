import re
from collections import Counter

from bgzip import BGZipWriter

from snpdb.models import VCF, Zygosity, Sample
from snpdb.vcf_export_utils import get_vcf_header_from_contigs, get_vcf_header_lines


def key_data_func(key, data):
    return data[key]


VARIANT_ID = 'variant_id'
VARIANT_PATH = 'variant_path'
VARIANT_GRID_INFO_DICT = {
    VARIANT_ID: {
        'type': 'Integer',
        'description': 'VariantGrid primary column',
        VARIANT_PATH: 'id',
        'key_data_func': key_data_func},
    # INFO fields for CNV
    "END": {
        'type': 'Integer',
        'description': 'Stop position of the interval',
    },
    "SVTYPE": {
        'type': 'String',
        'description': 'Type of structural variant',
    }
}


def qs_info_dict_field_values(qs, info_dict):
    args = {"locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq"}
    for data in info_dict.values():
        variant_path = data.get(VARIANT_PATH)
        if variant_path:
            args.add(variant_path)
    args = tuple(args)
    return qs.values(*args)


def write_qs_to_vcf_file_sort_alphabetically(qs, f, info_dict=None):
    if info_dict is None:
        info_dict = VARIANT_GRID_INFO_DICT
    header_lines = get_vcf_header_lines(info_dict=info_dict)

    sorted_query = qs.order_by("locus__contig__name", "locus__position")
    sorted_values = qs_info_dict_field_values(sorted_query, info_dict)
    return _write_sorted_values_to_vcf_file(header_lines, sorted_values, f, info_dict=info_dict)


def _write_sorted_values_to_vcf_file(header_lines, sorted_values, f, info_dict):
    symbolic_pattern = re.compile("<(.*)>")

    for line in header_lines:
        line_bytes = (line + '\n').encode()
        f.write(line_bytes)

    i = 0
    for data in sorted_values:
        chrom = data["locus__contig__name"]
        pos = data["locus__position"]
        ref = data["locus__ref__seq"]
        alt = data["alt__seq"]
        end = data["end"]
        info = {}

        if m := symbolic_pattern.match(alt):
            info["END"] = end
            info["SVTYPE"] = m.group(1)

        if info_dict:
            for info_name, info_data in info_dict.items():
                if variant_path := info_data.get("variant_path"):
                    func = info_data.get("key_data_func", key_data_func)
                    value = func(variant_path, data)
                    info[info_name] = value

        if info:
            info_str = ';'.join([f"{k}={v}" for k, v in info.items()])
        else:
            info_str = '.'
        line = '\t'.join(map(str, (chrom, pos, '.', ref, alt or ref, '.', '.', info_str)))
        line_bytes = (line + '\n').encode()
        f.write(line_bytes)
        i += 1

    return i


def write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f, info_dict):
    header_lines = get_vcf_header_from_contigs(genome_build, info_dict=info_dict)
    return _write_sorted_values_to_vcf_file(header_lines, sorted_values, f, info_dict=info_dict)


def vcf_export_to_file(vcf: VCF, exported_vcf_filename, original_qs=None, sample_name_func=None) -> dict[Sample, Counter]:
    """ Returns dict of zygosity counts written to file """
    if sample_name_func is None:
        def sample_name_func(s):
            return s.vcf_sample_name

    qs = vcf.get_variant_qs(original_qs)
    ca = vcf.cohort.cohort_genotype_collection.cohortgenotype_alias
    # Restrict to just this build (was returning multiple results due to GRCh37/hg19)
    qs = qs.filter(locus__contig__genomebuildcontig__genome_build=vcf.genome_build,
                   **{f"{ca}__filters__isnull": True})  # Somalier only uses PASS by default
    columns = ["id", "locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq",
               f"{ca}__samples_zygosity", f"{ca}__samples_allele_depth",
               f"{ca}__samples_read_depth", f"{ca}__samples_allele_frequency"]
    qs = qs.order_by("locus__contig__genomebuildcontig__order", "locus__position")

    vcf_format = "GT:AD:DP"
    samples = list(vcf.sample_set.order_by("pk"))
    sample_whitelist = [not s.no_dna_control for s in samples]  # Skip no DNA controls
    vcf_sample_names = [sample_name_func(s) for s, w in zip(samples, sample_whitelist) if w]
    header_lines = get_vcf_header_from_contigs(vcf.genome_build, samples=vcf_sample_names)
    sample_zygosity_count = [Counter() for _ in samples]
    empty = [None] * len(samples)

    with open(exported_vcf_filename, "wb") as raw:
        with BGZipWriter(raw) as f:
            for line in header_lines:
                line_bytes = (line + '\n').encode()
                f.write(line_bytes)

            values = qs.values_list(*columns)
            for pk, chrom, position, ref, alt, samples_zygosity, allele_depth, read_depth, allele_frequency in values:
                if allele_depth is None:
                    allele_depth = empty
                if read_depth is None:
                    read_depth = empty
                if allele_frequency is None:
                    allele_frequency = empty

                samples_list = []
                for i, (z, ad, dp, af) in enumerate(zip(samples_zygosity, allele_depth, read_depth, allele_frequency)):
                    if sample_whitelist[i]:
                        sample_zygosity_count[i][z] += 1
                        if z == Zygosity.UNKNOWN_ZYGOSITY:
                            sample = "./."
                        else:
                            gt = Zygosity.get_genotype(z)
                            if ad is not None:
                                if dp is not None and af is not None:
                                    if vcf.allele_frequency_percent:
                                        ref_depth = round(dp * (100 - af) / 100)
                                    else:
                                        ref_depth = round(dp * (1 - af))
                                else:
                                    ref_depth = '.'
                                ad = f"{ref_depth},{ad}"
                            else:
                                ad = '.'
                            if dp is None:
                                dp = '.'
                            sample = ":".join((str(s) for s in (gt, ad, dp)))
                        samples_list.append(sample)

                row = [chrom, str(position), str(pk), ref, alt or ref, '.', '.', '.', vcf_format] + samples_list
                line = '\t'.join(row)
                line_bytes = (line + '\n').encode()
                f.write(line_bytes)

    return {s: zc for s, zc, w in zip(samples, sample_zygosity_count, sample_whitelist) if w}
