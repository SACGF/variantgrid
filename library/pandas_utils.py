import os
import sys

import numpy as np
import pandas as pd


def get_columns_percent_dataframe(df, **kwargs):
    """ @param totals_column: (default = use sum of columns)
        @param percent_names: Rename names from 'col' => 'col %'
    Return a dataframe as a percentage of totals_column if provided, or sum of columns """

    totals_column = kwargs.get("totals_column")
    percent_names = kwargs.get("percent_names", True)

    percent_df = pd.DataFrame(index=df.index)
    columns = df.columns

    if totals_column:
        totals_series = df[totals_column]
        columns = columns - [totals_column]
    else:
        totals_series = df.sum(axis=1)

    for col in columns:
        new_col = col
        if percent_names:
            new_col = f"{new_col} %"
        multiplier = 100.0  # to get percent
        percent_df[new_col] = multiplier * df[col] / totals_series
    return percent_df


def get_rows_percent_dataframe(df):
    """ Return a dataframe as a percentage of sum of rows """
    row_sums = df.sum(axis=0)
    return 100.0 * df / row_sums


def get_total_percent_dataframe(df):
    """ Return a dataframe as a percentage of sum of rows """
    total = df.sum(axis=0).sum()
    return 100.0 * df / total


def df_handle_below_minimum_floats(df):

    def handle_if_below_min(series):
        if series.dtype == 'd':
            too_small_mask = abs(series) < sys.float_info.min
            series[too_small_mask] = sys.float_info.min

        return series

    return df.apply(handle_if_below_min, axis=0)


def nan_to_none(val):
    if np.isnan(val):
        val = None
    return val


def df_nan_to_none(df):
    return df.where((pd.notnull(df)), None)


def df_replace_nan(df, nan_replace=''):
    return df.where((pd.notnull(df)), nan_replace)


def read_csv_skip_header(fle, header='#', **kwargs):
    if os.stat(fle).st_size == 0:
        raise ValueError("File is empty")
    with open(fle) as f:
        pos = 0
        cur_line = f.readline()
        while cur_line.startswith(header):
            pos = f.tell()
            cur_line = f.readline()
        f.seek(pos)
        return pd.read_csv(f, **kwargs)
