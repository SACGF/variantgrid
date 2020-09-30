import operator
from functools import reduce

from library.date_utils import get_month_and_year, get_months_since, month_range
from seqauto.models import SequencingSample
import pandas as pd
import numpy as np


def get_sample_enrichment_kits_df():
    SEQUENCING_RUN_COL = "sample_sheet__sequencing_run"
    values_qs = SequencingSample.get_current().values(SEQUENCING_RUN_COL, "enrichment_kit__name")
    df = pd.DataFrame.from_records(values_qs)
    if SEQUENCING_RUN_COL in df.columns:
        sr = df[SEQUENCING_RUN_COL]

        year_month_series = pd.Series(index=df.index, dtype='i')
        year_series = pd.Series(index=df.index, dtype='i')

        for (i, val) in sr.items():
            run_date = val.split("_")[0]
            year_series[i] = run_date[:2]
            year_month_series[i] = run_date[:4]

        start_month, start_year = get_month_and_year(year_month_series.min())

        month_offset = pd.Series(index=df.index)
        for (i, year_month) in year_month_series.items():
            month, year = get_month_and_year(year_month)
            month_offset[i] = get_months_since(start_month, start_year, month, year)

        df.loc[:, "year"] = year_series
        df.loc[:, "year_month"] = year_month_series
        df.loc[:, "month_offset"] = month_offset
    return df


def year_formatter_start_to_end(start, end, _year_month_start):
    return [f"{i}" for i in range(start, end + 1)]


def year_month_formatter_start_to_end(start, end, year_month_start):
    start_month, start_year = get_month_and_year(year_month_start)
    return month_range(start_month, start_year, start, end)


def group_enrichment_kits_df(df, by_column, max_groups=None):
    """ returns (array of (enrichment_kit_name, data), labels)
        max_groups=10 gives 9 groups with everything else as "other" """
    LABELS_FOR_COLUMNS = {"year": year_formatter_start_to_end,
                          "month_offset": year_month_formatter_start_to_end}

    enrichment_kit_data = []
    labels = []

    if not df.empty:
        start = int(df[by_column].min())
        end = int(df[by_column].max())
        year_month_start = df["year_month"].min()

        array_size = int(end - start + 1)

        get_labels_from_start_to_end = LABELS_FOR_COLUMNS.get(by_column)
        if not get_labels_from_start_to_end:
            msg = "by_column must be one of %s" % ', '.join(LABELS_FOR_COLUMNS)
            raise ValueError(msg)
        labels = get_labels_from_start_to_end(start, end, year_month_start)

        for (enrichment_kit_name, enrichment_kit_df) in df.groupby("enrichment_kit__name"):
            array = [0] * array_size

            for value in enrichment_kit_df[by_column]:
                offset = int(value - start)
                array[offset] += 1

            enrichment_kit_data.append((enrichment_kit_name, array))

    if max_groups is not None and len(enrichment_kit_data) > max_groups:
        named_groups = max_groups - 1
        enrichment_kit_data_sum = [(name, array, sum(array)) for (name, array) in enrichment_kit_data]
        enrichment_kit_data_sum = list(sorted(enrichment_kit_data_sum, key=operator.itemgetter(2), reverse=True))
        enrichment_kit_data = []
        for (name, array, _) in enrichment_kit_data_sum[:named_groups]:
            enrichment_kit_data.append((name, array))

        other_sum = reduce(operator.add, [np.array(array) for (_, array, _) in enrichment_kit_data_sum[named_groups:]])
        enrichment_kit_data.append(("other", other_sum.tolist()))

    return enrichment_kit_data, labels
