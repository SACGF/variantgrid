import csv
import io
import math
import re
import string
from typing import Collection, Any, Optional, Callable


def pretty_label(label: str) -> str:
    label = label.replace('_', ' ')
    tidied = ''
    last_space = True
    for char in label:
        if last_space:
            char = char.upper()
            last_space = False
        if char == ' ':
            last_space = True
        tidied = tidied + char
    return tidied


def pretty_collection(collection: Collection[Any], to_string: Optional[Callable] = None) -> str:
    try:
        collection = sorted(collection)
    except:
        pass
    if to_string:
        collection = (to_string(item) for item in collection)

    return ", ".join(f'{item}' for item in collection)


def none_to_blank_string(s: Optional[str]) -> str:
    return s or ''


# don't think is still being used, would be passed into formatters
def upper(text: str) -> str:
    if text:
        text = str(text).upper()
    return text


def single_quote(s: Any) -> str:
    return f"'{s}'"


def double_quote(s: Any) -> str:
    return f'"{s}"'

def format_percent(number, is_unit=False) -> str:
    if is_unit:
        number *= 100
    return f"{format_significant_digits(number)}%"


trailing_zeros_strip = re.compile("(.*?[.][0-9]*?)(0+)$")


def format_significant_digits(a_number, sig_digits=3) -> str:
    if a_number == 0:
        return "0"
    rounded_number = round(a_number, sig_digits - int(math.floor(math.log10(abs(a_number)))) - 1)
    rounded_number_str = "{:.12f}".format(rounded_number)
    if match := trailing_zeros_strip.match(rounded_number_str):
        rounded_number_str = match.group(1)
        if rounded_number_str[-1] == '.':
            rounded_number_str = rounded_number_str[:-1]

    return rounded_number_str


def delimited_row(data: list, delimiter: str = ',') -> str:
    out = io.StringIO()
    writer = csv.writer(out, delimiter=delimiter)
    writer.writerow(data)
    return out.getvalue()


def clean_string(input_string: str) -> str:
    """ Removes non-printable characters, strips whitespace """
    return re.sub(f'[^{re.escape(string.printable)}]', '', input_string.strip())
