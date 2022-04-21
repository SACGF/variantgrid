import math
import re
from typing import Any

# any value other than weeks_gestation or months is treated as years
NON_FUZZY_AGE_RE = re.compile('([0-9]+)(weeks_gestation|months)?')


def patch_fuzzy_age(age: Any) -> str:
    """
    Converts ages to one of the following
    prenatal
    90+ (years)
    X0-X9 (years) - if above 0 aka prenatal and below 90

    Important - works on the combined age field, if you have historic data with age_units
    run patch_merge_age_units first
    :param age: The str to fuzzy up
    :return the fuzzy age, or just the input if it wasn't able to be processed
    """
    if age:
        age = str(age)
        if m := NON_FUZZY_AGE_RE.match(age):
            unit = m[2] or 'years'
            number_part = int(m[1])
            if unit == 'weeks_gestation':
                return 'prenatal'
            if unit == 'months':
                number_part = number_part / 12
                unit = 'years'

            if unit != 'years':
                raise ValueError(f'Unexpected unit {unit} in age fuzzier')

            if number_part >= 80:
                return '80+'
            else:
                base = math.floor(number_part / 10)
                if base == 0:
                    return '0-9'
                else:
                    return f'{base}0-{base}9'
    return age
