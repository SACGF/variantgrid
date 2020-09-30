import pandas as pd


def load_gene_coverage_df(gene_coverage_file):
    GENE_EXON_COLUMN = "Gene/Exon"
    UNDER_10X_COLUMN = "% bases <10x"
    UNDER_20X_COLUMN = "% bases <20x"
    OVER_100X_COLUMN = "% bases >100x"
    RENAME_COLUMS = {"Min Coverage": "min",
                     "Mean Coverage": "mean",
                     "Standard Deviation Coverage": "std_dev",
                     "% bases >20x": "percent_20x",
                     "% bases @ 0x": "percent_0x",
                     "Estimated Sensitivity of Detection": "sensitivity"}

    df = pd.read_csv(gene_coverage_file, sep='\t')
    # Split up into separate columns
    gene_exon_series = df[GENE_EXON_COLUMN].str.split(":")
    df["gene"] = gene_exon_series.str[0]
    df["transcript"] = gene_exon_series.str[1]

    if OVER_100X_COLUMN in df.columns:
        RENAME_COLUMS[OVER_100X_COLUMN] = "percent_100x"

    df = df.rename(columns=RENAME_COLUMS)

    # Coverage files currently generated inconsistency, sometimes using OVER
    # and sometimes under. Change them all to be over
    if UNDER_10X_COLUMN in df.columns:
        under_10x = df[UNDER_10X_COLUMN]
        df["percent_10x"] = 100.0 - under_10x
        del df[UNDER_10X_COLUMN]

    if UNDER_20X_COLUMN in df.columns:
        under_20x = df[UNDER_20X_COLUMN]
        df["percent_20x"] = 100.0 - under_20x
        del df[UNDER_20X_COLUMN]

    del df[GENE_EXON_COLUMN]
    return df
