from functools import partial

from library.utils.text_utils import format_significant_digits


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
    """ Shows falsey values (eg 0.0) or '.' as blank. Significant digits rather than '%g' - rare
        allele frequencies are small enough that '%g' switches to scientific notation """
    display_value = ""
    if val and val != missing_value:
        display_value = f"{format_significant_digits(val)}%"
    return display_value


def _get_formatters(source_in_percent, dest_in_percent, missing_value=None) -> list:
    formatters = []

    if dest_in_percent:
        if not source_in_percent:
            formatters.append(partial(convert_from_unit_to_percent, missing_value=missing_value))
        # Always run through format percent to add "%" etc
        formatters.append(partial(server_side_format_percent, missing_value=missing_value))
    else:
        if source_in_percent:
            formatters.append(partial(convert_from_percent_to_unit, missing_value=missing_value))
        # Do we want a unit formatter, e.g. to limit sig figures?
    return formatters


def format_af(value, source_in_percent, dest_in_percent, missing_value=None):
    formatters = _get_formatters(source_in_percent, dest_in_percent, missing_value=missing_value)
    for f in formatters:
        value = f(value)
    return value


def get_allele_frequency_formatter(source_in_percent, dest_in_percent, get_data_func=None, missing_value=None):
    """ A grid column renderer (@see snpdb.views.datatable_view.RichColumn) - get_data_func pulls the
        raw value out of the row where it isn't the column's own key (packed genotype columns) """
    formatters = _get_formatters(source_in_percent, dest_in_percent, missing_value=missing_value)

    def format_cell(cell):
        val = get_data_func(cell) if get_data_func else cell.value
        for f in formatters:
            val = f(val)
        return val

    return format_cell
