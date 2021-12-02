from functools import partial
from typing import List


def convert_from_percent_to_unit(percent, missing_value=None):
    if percent != missing_value:
        unit = percent / 100.0
    else:
        unit = missing_value
    return unit


def convert_from_unit_to_percent(unit, missing_value=None):
    if unit != missing_value:
        percent = unit * 100.0
    else:
        percent = missing_value
    return percent


def server_side_format_percent(val, missing_value=None):
    """ Shows falsey values (eg 0.0) or '.' as blank """
    display_value = ""
    if val and val != missing_value:
        display_value = f"{val:.3g}%"
    return display_value


def _get_formatters(source_in_percent, dest_in_percent, missing_value=None) -> List:
    formatters = []

    if dest_in_percent:
        if not source_in_percent:
            formatters.append(partial(convert_from_unit_to_percent, missing_value=missing_value))
        # Always run through format percent to add "%" etc
        formatters.append(partial(server_side_format_percent, missing_value=missing_value))
    else:
        if source_in_percent:
            formatters.append(partial(convert_from_percent_to_unit, missing_value=missing_value))
        # Do we want a unit formatter, eg to limit sig figures?
    return formatters


def format_af(value, source_in_percent, dest_in_percent, missing_value=None):
    formatters = _get_formatters(source_in_percent, dest_in_percent, missing_value=missing_value)
    for f in formatters:
        value = f(value)
    return value


def get_allele_frequency_formatter(source_in_percent, dest_in_percent, get_data_func=None, missing_value=None):
    formatters = _get_formatters(source_in_percent, dest_in_percent, missing_value=missing_value)

    def format_field(row, field):
        if get_data_func:
            val = get_data_func(row, field)
        else:
            val = row[field]
        for f in formatters:
            val = f(val)
        return val

    return format_field
