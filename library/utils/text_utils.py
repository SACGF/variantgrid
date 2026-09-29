import csv
import io
import math
import re
import string
from collections.abc import Callable, Collection
from typing import Any, Optional


def pretty_label(label: str) -> str:
    label = label.replace('_', ' ')
    tidied = ''
    last_space = True
    for char in label:
        if last_space:
            char = char.upper()
            last_space = False
        if char == ' ' or char == '-':
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
        parts = v.split(sep)
        if len(parts) != num_records:
            raise ValueError(f"Key {k!r} split into {len(parts)} parts, expected {num_records}")
        for i, v_part in enumerate(parts):
            dict_list[i][k] = v_part
    return dict_list


def limit_str(text: str, limit: int) -> str:
    if len(text) > limit:
        text = text[:limit] + "..."
    return text


def pretty_collection(collection: Collection[Any], to_string: Optional[Callable] = None) -> str:
    try:
        collection = sorted(collection)
    except Exception:
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
    rounded_number_str = f"{rounded_number:.12f}"
    if match := trailing_zeros_strip.match(rounded_number_str):
        rounded_number_str = match.group(1)
        if rounded_number_str[-1] == '.':
            rounded_number_str = rounded_number_str[:-1]

    return rounded_number_str


def delimited_row(data: list, delimiter: str = ',', include_new_line=True, **kwargs) -> str:
    # https://docs.python.org/3/library/csv.html#csv.writer
    # If csvfile is a file object, it should be opened with newline=''
    out = io.StringIO(newline='')
    if delimiter == '\t':
        writer = csv.writer(out, delimiter=delimiter, **kwargs)
    else:
        writer = csv.writer(out, delimiter=delimiter, **kwargs)
    writer.writerow(data)
    text = out.getvalue()
    if not include_new_line:
        text = text.rstrip()
    return text


def clean_string(input_string: str) -> str:
    if input_string is None:
        return ""
    """ Removes non-printable characters, strips whitespace """
    return re.sub(f'[^{re.escape(string.printable)}]', '', input_string.strip())


# Slack emoji codes used in notifications, health checks and settings (SLACK "emoji"). Unknown codes are left as is
_SLACK_EMOJI = {
    "arrow_down": "\u2b07\ufe0f",
    "bangbang": "\u203c\ufe0f",
    "blue_book": "\U0001f4d8",
    "cop": "\U0001f46e",
    "cowboy_hat_face": "\U0001f920",
    "crown": "\U0001f451",
    "cry": "\U0001f622",
    "currency_exchange": "\U0001f4b1",
    "dizzy_face": "\U0001f635",
    "dna": "\U0001f9ec",
    "email": "\u2709\ufe0f",
    "exploding_head": "\U0001f92f",
    "face_with_cowboy_hat": "\U0001f920",
    "female-doctor": "\U0001f469\u200d\u2695\ufe0f",
    "file_folder": "\U0001f4c1",
    "fire": "\U0001f525",
    "flags": "\U0001f38f",
    "floppy_disk": "\U0001f4be",
    "ghost": "\U0001f47b",
    "golfer": "\U0001f3cc\ufe0f",
    "green_book": "\U0001f4d7",
    "handshake": "\U0001f91d",
    "hospital": "\U0001f3e5",
    "hourglass_flowing_sand": "\u23f3",
    "male-doctor": "\U0001f468\u200d\u2695\ufe0f",
    "man_health_worker": "\U0001f468\u200d\u2695\ufe0f",
    "mouse": "\U0001f42d",
    "nerd_face": "\U0001f913",
    "neutral_face": "\U0001f610",
    "no_good": "\U0001f645",
    "open_file_folder": "\U0001f4c2",
    "orange_book": "\U0001f4d9",
    "package": "\U0001f4e6",
    "rage": "\U0001f621",
    "simple_smile": "\U0001f604",
    "skunk": "\U0001f9a8",
    "smile": "\U0001f604",
    "test_tube": "\U0001f9ea",
    "triangular_ruler": "\U0001f4d0",
    "tv": "\U0001f4fa",
    "warning": "\u26a0\ufe0f",
    "woman_health_worker": "\U0001f469\u200d\u2695\ufe0f",
}


def emoji_to_unicode(text_with_emojis: str) -> str:
    return re.sub(r":([a-z0-9_+-]+):", lambda m: _SLACK_EMOJI.get(m.group(1), m.group(0)), text_with_emojis)
