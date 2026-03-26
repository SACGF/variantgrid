import uuid

import pandas as pd

from patients.models_enums import Sex
from pedigree.ped.ped_file_utils import PED_COLUMNS


def write_ped_from_df(filename, df: pd.DataFrame, family_code=None):
    if family_code is None:
        family_code = str(uuid.uuid4())
    df["family"] = family_code
    df = df.fillna("0")
    df.to_csv(filename, sep='\t', index=False, header=False)


def write_trio_ped(filename, proband, proband_sex, father, father_affected, mother, mother_affected, family_code=None):
    # PED affection: 2=affected, 1=unaffected (standard numeric encoding)
    # PED sex: 1=male, 2=female
    sex_lookup = {Sex.MALE: "1", Sex.FEMALE: "2"}
    proband_sex = sex_lookup.get(proband_sex, "0")

    trio_data = [
        {'sample': proband, 'father': father, 'mother': mother, 'sex': proband_sex, 'affection': "2"},
        {'sample': father, 'sex': "1", 'affection': "2" if father_affected else "1"},
        {'sample': mother, 'sex': "2", 'affection': "2" if mother_affected else "1"},
    ]
    df = pd.DataFrame.from_records(trio_data, columns=PED_COLUMNS)
    write_ped_from_df(filename, df, family_code=family_code)


def write_unrelated_ped(filename, samples_names: list[str], family_code=None):
    unrelated_data = [{"sample": s} for s in samples_names]
    df = pd.DataFrame.from_records(unrelated_data, columns=PED_COLUMNS)
    write_ped_from_df(filename, df, family_code=family_code)
