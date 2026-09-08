import io
from collections import Counter
from typing import Optional

from bgzip import BGZipWriter

from library.genomics.vcf_enums import GeneLevelSymbolicAlt
from library.genomics.vcf_writer import VCFWriter, symbolic_alt_info
from snpdb.models import VCF, CohortGenotype, Sample, Zygosity
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
    # INFO fields for Symbolic alts
    "END": {
        'type': 'Integer',
        'description': 'Stop position of the interval',
    },
    "SVLEN": {
        'type': 'Integer',
        'description': 'Difference in length between REF and ALT alleles',
    },
    "SVTYPE": {
        'type': 'String',
        'description': 'Type of structural variant',
    }
}


def qs_info_dict_field_values(qs, info_dict):
    args = {"locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq", "end", "svlen"}
    for data in info_dict.values():
        variant_path = data.get(VARIANT_PATH)
        if variant_path:
            args.add(variant_path)
    args = tuple(args)
    return qs.values(*args)


def write_qs_to_vcf_file_sort_alphabetically(qs, f, info_dict=None) -> int:
    if info_dict is None:
        info_dict = VARIANT_GRID_INFO_DICT
    header_lines = get_vcf_header_lines(info_dict=info_dict)

    sorted_query = qs.order_by("locus__contig__name", "locus__position")
    sorted_values = qs_info_dict_field_values(sorted_query, info_dict)
    return _write_sorted_values_to_vcf_file(header_lines, sorted_values, f, info_dict=info_dict)


def _write_sorted_values_to_vcf_file(header_lines, sorted_values, f, info_dict, use_accession=False,
                                     samples=None) -> int:
    """
    :param samples: sample names written into every record as a dummy heterozygous call - for tools that
                    reject sites-only VCFs. Must match the samples the header_lines were built with
    :return: number of lines written
    """
    fmt = None
    sample_calls = None
    if samples:
        fmt = "GT"
        sample_calls = ["0/1"] * len(samples)

    if use_accession:
        chrom_key = "locus__contig__refseq_accession"
    else:
        chrom_key = "locus__contig__name"

    writer = VCFWriter(f, header_lines)
    i = 0
    for data in sorted_values:
        chrom = data[chrom_key]
        pos = data["locus__position"]
        ref = data["locus__ref__seq"]
        alt = data["alt__seq"]
        # A gene-level variant has a gene id where a coordinate goes, so a VCF record of it would be a
        # lie - and the one that reaches VEP is the expensive way to find that out. Callers keep them
        # out by contig (@see snpdb.gene_level_variants); this is what says so if one slips through.
        if alt.startswith("<") and GeneLevelSymbolicAlt.parse(alt):
            raise ValueError(f"Gene-level variant '{chrom}:{pos} {ref}>{alt}' cannot be written to a VCF")
        # END/SVLEN/SVTYPE are only needed (and only selected) for symbolic alts
        info = symbolic_alt_info(alt, svlen=data.get("svlen"), end=data.get("end"))

        if info_dict:
            for info_name, info_data in info_dict.items():
                if variant_path := info_data.get("variant_path"):
                    func = info_data.get("key_data_func", key_data_func)
                    value = func(variant_path, data)
                    info[info_name] = value

        writer.write_record(chrom, pos, ref, alt or ref, info=info or None,
                            fmt=fmt, sample_calls=sample_calls)
        i += 1

    return i


def write_contig_sorted_values_to_vcf_file(genome_build, sorted_values, f, info_dict, use_accession=False,
                                           samples=None) -> int:
    header_lines = get_vcf_header_from_contigs(genome_build, info_dict=info_dict, use_accession=use_accession,
                                               samples=samples)
    return _write_sorted_values_to_vcf_file(header_lines, sorted_values, f, info_dict=info_dict,
                                            use_accession=use_accession, samples=samples)


SOMALIER_GT_FORMAT = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
SOMALIER_AD_FORMAT = '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">'


def _missing_to_none(value):
    """ CohortGenotype stores 'no value' as -1, so 'is not None' isn't enough """
    if value is None or value == CohortGenotype.MISSING_NUMBER_VALUE:
        return None
    return value


SOMALIER_FLIPPED_GENOTYPE = {"0/0": "1/1", "1/1": "0/0"}


def _somalier_alleles_flipped(ref: str, alt: str) -> bool:
    """ somalier keeps each site's two alleles in alphabetical order (A = min, B = max) and reads GT
        and AD against that pair, not against the record's REF/ALT - so a record whose ALT sorts
        first is read inside out. Measured on v0.2.12 and v0.3.4: 200 hom-ref calls at real sites
        come back as 133 hom-alt. Emitting the alleles in the site's own order is what makes somalier
        agree with the genotypes we imported. """
    return alt < ref


def _allele_depths(vcf: VCF, alt_depth, read_depth, allele_frequency) -> Optional[tuple[int, int]]:
    """ (ref, alt) depths. We only store the alt depth, so ref comes from DP or AF; None when neither
        gets us there - the caller writes '.', which somalier reads as 0,0 and its own --min-depth
        drops for that sample """
    alt_depth = _missing_to_none(alt_depth)
    if alt_depth is None:
        return None
    read_depth = _missing_to_none(read_depth)
    if read_depth is not None:
        return max(0, read_depth - alt_depth), alt_depth

    allele_frequency = _missing_to_none(allele_frequency)
    if not allele_frequency:
        return None
    if vcf.allele_frequency_percent:
        allele_frequency /= 100
    return max(0, round(alt_depth * (1 - allele_frequency) / allele_frequency)), alt_depth


def vcf_export_to_file(vcf: VCF, exported_vcf_filename, original_qs=None, sample_name_func=None) -> dict[Sample, Counter]:
    """ Writes a VCF for 'somalier extract'. Returns dict of zygosity counts written to file.

        Somalier decides how to genotype from the header (v0.2.12 get_ref_alt_counts): a FORMAT AD line
        means it re-genotypes every sample from the depths, and applies its own QC at relate time
        (--min-depth 7, --min-ab 0.3); with no AD line it trusts GT at pseudo-depths. VariantGrid takes
        the incoming GT as called, however low the depth, so we hand over real depths when the import
        recorded them and let somalier do the QC - and GT alone when it didn't, since declaring AD we
        can't fill in would zero out every sample. Every record goes out in the site's own allele order
        (@see _somalier_alleles_flipped). """
    if sample_name_func is None:
        def sample_name_func(s):
            return s.vcf_sample_name

    # A ref depth is derivable from DP or AF; without one of them there's nothing honest to put in AD
    write_allele_depth = bool(vcf.allele_depth_field and (vcf.read_depth_field or vcf.allele_frequency_field))

    qs = vcf.get_variant_qs(original_qs)
    ca = vcf.cohort.cohort_genotype_collection.cohortgenotype_alias
    # Restrict to just this build (was returning multiple results due to GRCh37/hg19)
    qs = qs.filter(locus__contig__genomebuildcontig__genome_build=vcf.genome_build,
                   **{f"{ca}__filters__isnull": True})  # Somalier only uses PASS by default
    columns = ["id", "locus__contig__name", "locus__position", "locus__ref__seq", "alt__seq",
               f"{ca}__samples_zygosity"]
    if write_allele_depth:
        columns += [f"{ca}__samples_allele_depth", f"{ca}__samples_read_depth",
                    f"{ca}__samples_allele_frequency"]
    qs = qs.order_by("locus__contig__genomebuildcontig__order", "locus__position")

    if write_allele_depth:
        vcf_format = "GT:AD"
        formats = [SOMALIER_GT_FORMAT, SOMALIER_AD_FORMAT]
        unknown_call = "./.:."
    else:
        vcf_format = "GT"
        formats = [SOMALIER_GT_FORMAT]
        unknown_call = "./."

    samples = list(vcf.sample_set.order_by("pk"))
    sample_whitelist = [not s.no_dna_control for s in samples]  # Skip no DNA controls
    vcf_sample_names = [sample_name_func(s) for s, w in zip(samples, sample_whitelist) if w]
    # use_accession=False so the ##contig IDs match the CHROM we write and the 'nochr' sites file
    header_lines = get_vcf_header_from_contigs(vcf.genome_build, samples=vcf_sample_names,
                                               use_accession=False, formats=formats)
    sample_zygosity_count = [Counter() for _ in samples]
    empty = [None] * len(samples)

    with open(exported_vcf_filename, "wb") as raw:
        with BGZipWriter(raw) as bgzip_f:
            # bgzip is binary; wrap as text so VCFWriter only ever deals with str
            f = io.TextIOWrapper(bgzip_f, encoding="utf-8", write_through=True)
            writer = VCFWriter(f, header_lines)

            for row in qs.values_list(*columns):
                pk, chrom, position, ref, alt, samples_zygosity = row[:6]
                if write_allele_depth:
                    allele_depth, read_depth, allele_frequency = (v if v is not None else empty for v in row[6:])
                else:
                    allele_depth = read_depth = allele_frequency = empty

                alt = alt or ref
                flipped = _somalier_alleles_flipped(ref, alt)

                samples_list = []
                for i, (z, ad, dp, af) in enumerate(zip(samples_zygosity, allele_depth, read_depth, allele_frequency)):
                    if sample_whitelist[i]:
                        sample_zygosity_count[i][z] += 1
                        if z == Zygosity.UNKNOWN_ZYGOSITY:
                            sample = unknown_call
                        else:
                            sample = Zygosity.get_genotype(z)
                            if flipped:
                                sample = SOMALIER_FLIPPED_GENOTYPE.get(sample, sample)
                            if write_allele_depth:
                                depths = _allele_depths(vcf, ad, dp, af)
                                if depths is None:
                                    sample += ":."
                                else:
                                    ref_depth, alt_depth = depths
                                    if flipped:
                                        ref_depth, alt_depth = alt_depth, ref_depth
                                    sample += f":{ref_depth},{alt_depth}"
                        samples_list.append(sample)

                if flipped:
                    ref, alt = alt, ref
                writer.write_record(chrom, position, ref, alt, vcf_id=pk,
                                    fmt=vcf_format, sample_calls=samples_list)

            # flush + detach so the BGZipWriter (closed by its 'with') is closed exactly once
            f.flush()
            f.detach()

    return {s: zc for s, zc, w in zip(samples, sample_zygosity_count, sample_whitelist) if w}
