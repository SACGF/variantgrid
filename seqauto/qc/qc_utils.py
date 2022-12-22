import os
import re

from library.utils.file_utils import name_from_filename


def meta_data_file(filename, meta_data_file_format):
    basedir = os.path.dirname(filename)
    name = name_from_filename(filename)
    meta_data_relative_filename = meta_data_file_format % name
    return os.path.join(basedir, meta_data_relative_filename)


def match_patterns_in_file(source, patterns, mandatory):
    data = {}
    i = 0
    for line in source:
        i += 1
        for name, pattern in patterns.items():
            match_obj = re.search(pattern, line)
            if match_obj:
                data[name] = match_obj
                break

    if mandatory:
        num_lines = i
        for name, pattern in patterns.items():
            if not data.get(name):
                error_str = f"Pattern: {name} = {pattern} had no matches (searched {num_lines} lines)"
                raise ValueError(error_str)

    return data
