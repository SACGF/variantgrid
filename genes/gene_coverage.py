import pandas as pd


def load_gene_coverage_df(gene_coverage_file):
    GENE_TRANSCRIPT_VERSION_COLUMN = "Name"
    RENAME_COLUMNS = {
        'Size': "size",
        'Mean': "mean",
        'StdDev': "std_dev",
        'Max': "max",
        'Min': "min",
        '% >=1X': "percent_1x",
        '% >=2X': "percent_2x",
        '% >=5X': "percent_5x",
        '% >=10X': "percent_10x",
        '% >=15X': "percent_15x",
        '% >=20X': "percent_20x",
        '% >=25X': "percent_25x",
        '% >=30X': "percent_30x",
        '% >=40X': "percent_40x",
        '% >=50X': "percent_50x",
        '% >=60X': "percent_60x",
        '% >=80X': "percent_80x",
        '% >=100X': "percent_100x",
        '% >=150X': "percent_150x",
        '% >=200X': "percent_200x",
        '% >=250X': "percent_250x",
    }

    df = pd.read_csv(gene_coverage_file, sep='\t')
    # Split up into separate columns
    gene_transcript_version_series = df[GENE_TRANSCRIPT_VERSION_COLUMN].str.split(";")
    df["original_gene_symbol"] = gene_transcript_version_series.str[0]
    df["original_transcript"] = gene_transcript_version_series.str[1]

    del df[GENE_TRANSCRIPT_VERSION_COLUMN]
    df = df.rename(columns=RENAME_COLUMNS)
    return df
