import uuid
from typing import List

import pandas as pd

from library.utils import invert_dict
from pedigree.ped.ped_file_utils import SEX_LOOKUP, PED_COLUMNS


def write_ped_from_df(filename, df: pd.DataFrame, family_code=None):
    if family_code is None:
        family_code = str(uuid.uuid4())
    df["family"] = family_code
    df = df.fillna("0")
    df.to_csv(filename, sep='\t', index=False, header=False)


def write_trio_ped(filename, proband, proband_sex, father, father_affected, mother, mother_affected, family_code=None):
    sex_lookup = invert_dict(SEX_LOOKUP)
    proband_sex = sex_lookup.get(proband_sex, "0")

    trio_data = [
        {'sample': proband, 'father': father, 'mother': mother, 'sex': proband_sex, 'affection': "1"},
        {'sample': father, 'sex': "1", 'affection': "1" if father_affected else "0"},
        {'sample': mother, 'sex': "2", 'affection': "1" if mother_affected else "0"},
    ]
    df = pd.DataFrame.from_records(trio_data, columns=PED_COLUMNS)
    write_ped_from_df(filename, df, family_code=family_code)


def write_unrelated_ped(filename, samples_names: List[str], family_code=None):
    unrelated_data = [{"sample": s} for s in samples_names]
    df = pd.DataFrame.from_records(unrelated_data, columns=PED_COLUMNS)
    write_ped_from_df(filename, df, family_code=family_code)
