"""
Code by Mike T - http://stackoverflow.com/a/8150329

"""
import itertools
from struct import pack

import numpy as np


def prepare_binary(f, dat):
    """ f: file like object,
        dat: numpy array

        # @see https://www.postgresql.org/docs/10/static/datatype-numeric.html
        dtype = ([('variant_id', 'i4'), ('zygosity', 'S1'), ('allele_frequency', 'f4')])
    """

    pgcopy_dtype = [('num_fields', '>i2')]
    num_fields = len(dat.dtype)

    for field, dtype in dat.dtype.descr:
        pgcopy_dtype += [(field + '_length', '>i4'),
                         (field, dtype.replace('<', '>'))]
    pgcopy = np.empty(dat.shape, pgcopy_dtype)
    pgcopy['num_fields'] = num_fields

    for i in range(len(dat.dtype)):
        field = dat.dtype.names[i]
        pgcopy[field + '_length'] = dat.dtype[i].alignment
        pgcopy[field] = dat[field]

    f.write(pack('!11sii', b'PGCOPY\n\377\r\n\0', 0, 0))
    f.write(pgcopy.tostring())  # all rows
    f.write(pack('!h', -1))  # file trailer


def postgres_arrays(array):
    return "{%s}" % ','.join([str(s) if s is not None else "NULL" for s in array])


def split_postgres_arrays(pg_array_str, missing_fill_value=None):
    if pg_array_str is None:
        return itertools.cycle([missing_fill_value])
    return pg_array_str[1:-1].split(",")
