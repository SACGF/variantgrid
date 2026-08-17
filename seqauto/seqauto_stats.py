import logging
import operator
import re
from functools import reduce
from typing import Optional

import numpy as np
import pandas as pd
from django.conf import settings

from library.utils.date_utils import get_months_since, month_range, parse_yymm
from seqauto.models import SequencingRun, SequencingSample

SEQUENCING_RUN_DATE_PATTERN = re.compile(r"^(\d{4})\d{2}[_-]")


def get_sequencing_run_yymm(sequencing_run_name) -> Optional[str]:
    """ Illumina sequencing runs start with a YYMMDD date - returns the YYMM digits (None if unparseable) """
    if not isinstance(sequencing_run_name, str):
        return None
    sr_name = SequencingRun.get_original_illumina_sequencing_run(sequencing_run_name)
    if m := SEQUENCING_RUN_DATE_PATTERN.match(sr_name):
        return m.group(1)
    return None


def get_sample_enrichment_kits_df():
    if settings.SEQAUTO_FAKE_SAMPLE_ENRICHMENT_KITS_DF:
        logging.info("Loading FAKE sample enrichment kits DF: %s", settings.SEQAUTO_FAKE_SAMPLE_ENRICHMENT_KITS_DF)
        df = pd.read_csv(settings.SEQAUTO_FAKE_SAMPLE_ENRICHMENT_KITS_DF)
        return df

    SEQUENCING_RUN_COL = "sample_sheet__sequencing_run"
    values_qs = SequencingSample.get_current().values(SEQUENCING_RUN_COL, "enrichment_kit__name")
    df = pd.DataFrame.from_records(values_qs)
    # May be no sequencing runs - in which case skip
    if SEQUENCING_RUN_COL not in df.columns:
        return df

    years = {}
    year_months = {}
    for i, sequencing_run_name in df[SEQUENCING_RUN_COL].items():
        if yymm := get_sequencing_run_yymm(sequencing_run_name):
            years[i] = int(yymm[:2])
            year_months[i] = int(yymm)
        else:
            logging.warning("Skipping sequencing run %s - no YYMMDD date at start of name", sequencing_run_name)

    # Only keep rows we could date, so callers always see year/year_month/month_offset columns
    df = df.loc[list(year_months)].copy()
    if df.empty:
        return df

    year_month_series = pd.Series(year_months, dtype=int)
    start_month, start_year = parse_yymm(year_month_series.min())

    month_offsets = {}
    for i, year_month in year_month_series.items():
        month, year = parse_yymm(year_month)
        month_offsets[i] = get_months_since(start_month, start_year, month, year)

    df["year"] = pd.Series(years, dtype=int)
    df["year_month"] = year_month_series
    df["month_offset"] = pd.Series(month_offsets, dtype=int)
    return df


def year_formatter_start_to_end(start, end, _year_month_start):
    return [f"{i}" for i in range(start, end + 1)]


def year_month_formatter_start_to_end(start, end, year_month_start):
    start_month, start_year = parse_yymm(year_month_start)
    return month_range(start_month, start_year, start, end)


def group_enrichment_kits_df(df, by_column, max_groups=None, max_years=None):
    """ returns (array of (enrichment_kit_name, data), labels)
        max_groups=10 gives 9 groups with everything else as "other" """
    LABELS_FOR_COLUMNS = {"year": year_formatter_start_to_end,
                          "month_offset": year_month_formatter_start_to_end}

    enrichment_kit_data = []
    labels = []

    if not df.empty:
        if max_years is not None:
            last_year = df["year"].max()
            year_idx = df["year"] >= last_year - max_years
            df = df[year_idx]

        start = int(df[by_column].min())
        end = int(df[by_column].max())
        year_month_start = df["year_month"].min()

        array_size = int(end - start + 1)

        get_labels_from_start_to_end = LABELS_FOR_COLUMNS.get(by_column)
        if not get_labels_from_start_to_end:
            label_options = ', '.join(LABELS_FOR_COLUMNS)
            raise ValueError(f"by_column must be one of {label_options}")
        labels = get_labels_from_start_to_end(start, end, year_month_start)

        for enrichment_kit_name, enrichment_kit_df in df.groupby("enrichment_kit__name"):
            array = [0] * array_size

            for value in enrichment_kit_df[by_column]:
                offset = int(value - start)
                array[offset] += 1

            enrichment_kit_data.append((enrichment_kit_name, array))

    if max_groups is not None and len(enrichment_kit_data) > max_groups:
        named_groups = max_groups - 1
        enrichment_kit_data_sum = [(name, array, sum(array)) for name, array in enrichment_kit_data]
        enrichment_kit_data_sum = sorted(enrichment_kit_data_sum, key=operator.itemgetter(2), reverse=True)
        enrichment_kit_data = []
        for name, array, _ in enrichment_kit_data_sum[:named_groups]:
            enrichment_kit_data.append((name, array))

        other_sum = reduce(operator.add, [np.array(array) for _, array, _ in enrichment_kit_data_sum[named_groups:]])
        enrichment_kit_data.append(("other", other_sum.tolist()))

    return enrichment_kit_data, labels
