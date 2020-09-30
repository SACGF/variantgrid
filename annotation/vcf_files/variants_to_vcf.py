import os
import logging

from annotation.annotation_version_querysets import get_unannotated_variants_qs
from annotation.vcf_files.vcf_export_utils import get_vcf_header_from_contigs, get_vcf_header_lines
from library.file_utils import mk_path_for_file
from snpdb.models.models_genome import GenomeBuild


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


def write_unannotated_variants_to_vcf_file(genome_build, annotation_range_lock, f):
    kwargs = {}
    if annotation_range_lock:
        kwargs["min_variant_id"] = annotation_range_lock.min_variant.pk
        kwargs["max_variant_id"] = annotation_range_lock.max_variant.pk

    annotation_version = annotation_range_lock.version.get_any_annotation_version()
    qs = get_unannotated_variants_qs(annotation_version, **kwargs)
    qs = qs.order_by("locus__contig__genomebuildcontig__order", "locus__position")
    sorted_values = qs.values("id", "locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq")
    return write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f, info_dict=VARIANT_GRID_INFO_DICT)


def unannotated_variants_to_vcf(genome_build: GenomeBuild, vcf_filename, annotation_range_lock):
    logging.info("unannotated_variants_to_vcf()")
    if os.path.exists(vcf_filename):
        raise ValueError(f"Don't want to overwrite '{vcf_filename}' which already exists!")
    mk_path_for_file(vcf_filename)
    with open(vcf_filename, 'wb') as f:
        return write_unannotated_variants_to_vcf_file(genome_build, annotation_range_lock, f)
