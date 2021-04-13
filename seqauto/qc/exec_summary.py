"""
#995 - rewrite 28/10/2019 - read from CSV instead of human file
"""

import pandas as pd


def load_exec_summary(klass, exec_summary_filename):
    QCEXEC_FIELDS = {
        '%Duplication': 'percent_duplication',
        '%ErrorRate': 'percent_error_rate',
        '%Indels_dbSNP': 'indels_dbsnp_percent',
        '%ReadEnrich': 'percent_read_enrichment',
        '%Reads': 'percent_reads',
        '%SNPs_dbSNP': 'snp_dbsnp_percent',
        '%SoftClipRate': 'percent_softclip',
        '%UniformCov': 'uniformity_of_coverage',
        '%bp10X_GOI': 'percent_10x_goi',
        '%bp20X_GOI': 'percent_20x_goi',
        '%bp20x_Kit': 'percent_20x_kit',
        '%map2diff_chr': 'percent_map_to_diff_chr',
        'DedupReads': 'deduplicated_reads',
        'Indels': 'number_indels',
        'MeanDepth_GOI': 'mean_coverage_across_genes',
        'MeanDepth_Kit': 'mean_coverage_across_kit',
        'MedianInsert': 'median_insert',
        'Reads': 'reads',
        'SNPs': 'number_snps',
        'SampleID_LOD': 'sample_id_lod',
        'SexMatch': 'sex_match',
        'Ts/Tv': 'ts_to_tv_ratio',
    }
    type_converters = {f.name: f.get_internal_type() for f in klass._meta.fields}

    df = pd.read_csv(exec_summary_filename, sep='\t', index_col='Name')
    values = df["Value"]
    sample_name = values["Sample"]

    exec_data = {}
    for k, v in values.items():
        f = QCEXEC_FIELDS.get(k)
        if f and not pd.isna(v):
            internal_type = type_converters[f]
            # Integers come out as eg "42.0" so need to convert to float first then int
            if internal_type in ('FloatField', 'IntegerField'):
                v = float(v)
            if internal_type == 'IntegerField':
                v = int(v)
            exec_data[f] = v

    reference_range = {}  # No gold in new QC Summaries
    # for k, v in df[gold].items():
    #     if not pd.isna(v):
    #         if k == "MeanDepth_Kit":
    #             if m := re.match(r"> (\d+)x", v):
    #                 min_mean = m.group(1)
    #                 reference_range["min_mean_coverage_across_kit"] = float(min_mean)
    #         else:
    #             f = QCEXEC_FIELDS.get(k)
    #             if f:
    #                 (range_min, range_max) = v.split(" - ")
    #                 reference_range[f] = (float(range_min), float(range_max))

    data = {"sample": sample_name,
            "exec_data": exec_data,
            "reference_range": reference_range}
    return data
