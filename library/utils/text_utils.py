import csv
import io
import math
import re
import string
from typing import Collection, Any, Optional, Callable

from rich.text import Text


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
        tidied += char
    return tidied


def join_with_commas_and_ampersand(items: list[str], final_sep: str = "&") -> str:
    if len(items) == 0:
        return ""
    elif len(items) == 1:
        return items[0]
    elif len(items) == 2:
        return f" {final_sep} ".join(items)
    else:
        comma_sep = ", ".join(items[0:-1])
        return f" {final_sep} ".join([comma_sep, items[-1]])


def split_dict_multi_values(data: dict[str, str], sep: str) -> list[dict[str, str]]:
    any_value = next(iter(data.values()))
    num_records = len(any_value.split(sep))
    dict_list = [{} for _ in range(num_records)]
    for k, v in data.items():
        for i, v_part in enumerate(v.split(sep)):
            dict_list[i][k] = v_part
    return dict_list


def limit_str(text: str, limit: int) -> str:
    if len(text) > limit:
        text = text[:limit] + "..."
    return text


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


class VCFDialect(csv.Dialect):
    delimiter = ','
    lineterminator = "\n"
    quoting = csv.QUOTE_NONE


def delimited_row(data: list, delimiter: str = ',', dialect=None) -> str:
    kwargs = {}
    if dialect is not None:
        kwargs["dialect"] = dialect
    # https://docs.python.org/3/library/csv.html#csv.writer
    # If csvfile is a file object, it should be opened with newline=''
    out = io.StringIO(newline='')
    writer = csv.writer(out, delimiter=delimiter, **kwargs)
    writer.writerow(data)
    return out.getvalue()


def clean_string(input_string: str) -> str:
    if input_string is None:
        return ""
    """ Removes non-printable characters, strips whitespace """
    return re.sub(f'[^{re.escape(string.printable)}]', '', input_string.strip())


def emoji_to_unicode(text_with_emojis) -> str:
    # To see available emojis in Rich - python3 -m rich.emoji
    _replace = {
        ":male-doctor:": ":man_health_worker:",
        ":female-doctor:": ":woman_health_worker:",
        ":face_with_cowboy_hat:": ":cowboy_hat_face:",
        ":simple_smile:": ":smile:"
    }
    for old, new in _replace.items():
        text_with_emojis = text_with_emojis.replace(old, new)

    s = Text.from_markup(text_with_emojis)
    return str(s)
