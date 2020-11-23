from snpdb.vcf_export_utils import get_vcf_header_from_contigs, get_vcf_header_lines


def key_data_func(key, data):
    return data[key]


VARIANT_ID = 'variant_id'
VARIANT_PATH = 'variant_path'
VARIANT_GRID_INFO_DICT = {VARIANT_ID: {'number': 1,
                                       'type': 'Integer',
                                       'description': 'VariantGrid primary column',
                                       VARIANT_PATH: 'id',
                                       'key_data_func': key_data_func}}


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
    for line in header_lines:
        line_bytes = (line + '\n').encode()
        f.write(line_bytes)

    i = 0
    for data in sorted_values:
        chrom = data["locus__contig__name"]
        pos = data["locus__position"]
        ref = data["locus__ref__seq"]
        alt = data["alt__seq"]

        if info_dict:
            infos_list = []
            for info_name, info_data in info_dict.items():
                variant_path = info_data["variant_path"]
                func = info_data.get("key_data_func", key_data_func)
                value = func(variant_path, data)
                infos_list.append(f"{info_name}={value}")
            info = ';'.join(infos_list)
        else:
            info = '.'
        line = '\t'.join(map(str, (chrom, pos, '.', ref, alt or ref, '.', '.', info)))
        line_bytes = (line + '\n').encode()
        f.write(line_bytes)
        i += 1

    return i


def write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f, info_dict):
    header_lines = get_vcf_header_from_contigs(genome_build, info_dict=info_dict)
    return _write_sorted_values_to_vcf_file(header_lines, sorted_values, f, info_dict=info_dict)
