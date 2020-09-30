"""
#995 - rewrite 28/10/2019 - read from CSV instead of human file
"""

import re

import pandas as pd

EXEC_STATS_LOOKUP = [("% of bases with >20x coverage across kit", "percent_20x_kit"),
                     ("% bases >20x depth across kit", "percent_20x_kit"),
                     ("% of bases with >500x coverage across GOI", 'percent_500x'),
                     ("% bases >500x depth across GOI", 'percent_500x'),

                     ("% of bases with >250x coverage across GOI", 'percent_250x'),
                     ("% bases >250x depth across GOI", 'percent_250x'),

                     ("% of bases with >20x coverage across GOI", 'percent_20x'),
                     ("% bases >20x depth across GOI", 'percent_20x'),

                     ("% of bases with >10x coverage across GOI", 'percent_10x'),
                     ("% bases >10x depth across GOI", 'percent_10x'),

                     # Could be either of these 2
                     ("Mean coverage across GOI", 'mean_coverage_across_genes'),
                     ("Mean depth across GOI", 'mean_coverage_across_genes'),

                     # Could be either of these 3
                     ("Mean coverage across kit", 'mean_coverage_across_kit'),
                     ("Mean depth of coverage across kit", 'mean_coverage_across_kit'),
                     ("Mean depth across kit", 'mean_coverage_across_kit'),

                     ("Uniformity of Coverage", 'uniformity_of_coverage'),
                     ("Read enrichment (%)", 'percent_read_enrichment'),
                     ("Duplicated alignable reads", 'duplicated_alignable_reads'),
                     # https://bitbucket.org/sacgf/tau/issues/181/exec-summary-mean-insert-size-is-wrongly
                     # Old QC said MEAN but meant median, so we'll take both and store as median...
                     ("Mean insert", 'median_insert'),
                     ("Median insert", 'median_insert'),
                     ("Ts/Tv", 'ts_to_tv_ratio'),
                     ("Number SNPs", 'number_snps'),
                     ("%SNPs in dbSNP", 'snp_dbsnp_percent'),
                     ("Number indels", 'number_indels'),
                     ("%indels in dbSNP", 'indels_dbsnp_percent')]


def load_exec_summary(exec_summary_file):
    QCEXEC_FIELDS = {
        '%bp20x_GOI': 'percent_20x',
        '%bp10x_GOI': 'percent_10x',
        'MeanDepth_GOI': 'mean_coverage_across_genes',
        'MeanDepth_Kit': 'mean_coverage_across_kit',
        '%UniformCov': 'uniformity_of_coverage',
        '%ReadEnrich': 'percent_read_enrichment',
        '%Duplication': 'duplicated_alignable_reads',
        'MedianInsert': 'median_insert',
        'Ts/Tv': 'ts_to_tv_ratio',
        'SNPs': 'number_snps',
        '%SNPs_dbSNP': 'snp_dbsnp_percent',
        'Indels': 'number_indels',
        '%Indels_dbSNP': 'indels_dbsnp_percent'
    }

    df = pd.read_csv(exec_summary_file, sep='\t', index_col='Sample')
    assert len(df.columns) == 2, f"QCExec summary not a 2 column TSV, columns: {df.columns}"
    (sample_name, gold) = df.columns
    assert gold == "Gold", f"QCExec summary missing gold column in columns: {df.columns}"

    exec_data = {}
    for k, v in df[sample_name].items():
        f = QCEXEC_FIELDS.get(k)
        if f:
            exec_data[f] = float(v)

    reference_range = {}
    for k, v in df[gold].items():
        if not pd.isna(v):
            if k == "MeanDepth_Kit":
                if m := re.match(r"> (\d+)x", v):
                    min_mean = m.group(1)
                    reference_range["min_mean_coverage_across_kit"] = float(min_mean)
            else:
                f = QCEXEC_FIELDS.get(k)
                if f:
                    (range_min, range_max) = v.split(" - ")
                    reference_range[f] = (float(range_min), float(range_max))

    data = {"sample": sample_name,
            "exec_data": exec_data,
            "reference_range": reference_range}
    return data
