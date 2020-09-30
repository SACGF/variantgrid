import re

import numpy as np


# fastqc output from latest version has different numbering system
# creating this error - ValueError: invalid literal for int() with base 10: '100-101'
# will handle the error when presented and return none
def parse_sequence_length(lines):
    try:
        last_value = int(lines[-1].split()[0])
    except ValueError:
        return None
    else:
        array = np.zeros(last_value+1)
        for line in lines:
            if line.startswith("#"):
                continue
            (length, count) = line.split()
            array[int(length)] = count

    return array

def parse_basic_statistics(lines):
    basic_stats = {}

    type_converters = {
        "Total Sequences": int,
        "Sequences flagged as poor quality": int,
        "GC": int,
    }

    for line in lines:
        if line.startswith("#"):
            continue
        (key, value) = line.split('\t')
        key = key.replace("%", "")  # %GC => GC
        conversion_func = type_converters.get(key, str)
        basic_stats[key] = conversion_func(value)
    return basic_stats

#Section headers
BASIC_STATISTICS = "Basic Statistics"
SEQUENCE_LENGTH_DISTRIBUTION = "Sequence Length Distribution"

# Add lines here
SECTION_HANDLERS = {BASIC_STATISTICS: parse_basic_statistics,
                    SEQUENCE_LENGTH_DISTRIBUTION: parse_sequence_length}

def read_fastqc_data(fastqc_data_filename, custom_handlers=None):
    """ custom_handlers : dict {"section header" : function passed lines, returns parsed data}
        Returns a dict of section header : parsed data """

    handlers = SECTION_HANDLERS.copy()
    if custom_handlers:
        handlers.update(custom_handlers)

    section_pattern = re.compile(r">>(.*)\t")
    fastqc_data = {}

    with open(fastqc_data_filename) as f:
        lines = []
        section = None
        for line in f:
            line = line.rstrip()
            if line.startswith(">>"):
                if line[2:].startswith("END_MODULE"):
                    if section in handlers:
                        fastqc_data[section] = handlers[section](lines)
                else:
                    if m := section_pattern.match(line):
                        token = m.group(1)
                        section = token
                        lines = []
            else:
                lines.append(line)

    return fastqc_data
