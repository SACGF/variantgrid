import math
import re

from variantclassification.enums import SpecialEKeys
from variantclassification.models import PatchMeta


def patch_merge_age_units(patch_meta: PatchMeta):
    if patch_meta.intersection_modified({SpecialEKeys.AGE, SpecialEKeys.AGE_UNITS}):
        age = patch_meta.get(SpecialEKeys.AGE)
        age_units = patch_meta.get(SpecialEKeys.AGE_UNITS)

        combined_age = f'{age or ""}{age_units or ""}'
        if age != combined_age:
            patch_meta.patch_value(SpecialEKeys.AGE, combined_age)
            patch_meta.remove_patch_value(SpecialEKeys.AGE_UNITS)


NON_FUZZY_AGE_RE = re.compile('([0-9]+)(weeks_gestation|months)?')


def patch_fuzzy_age(patch_meta: PatchMeta):
    """
    Converts ages to one of the following
    prenatal
    90+ (years)
    X0-X9 (years) - if above 0 aka prenatal and below 90

    Important - works on the combined age field, if you have historic data with age_units
    run patch_merge_age_units first
    :param patch_meta: The patch to update
    """
    if patch_meta.is_modified(SpecialEKeys.AGE):
        age = str(patch_meta.get(SpecialEKeys.AGE))
        if age:
            if m := NON_FUZZY_AGE_RE.match(age):
                unit = m[2] or 'years'
                number_part = int(m[1])
                if unit == 'weeks_gestation':
                    patch_meta.patch_value(SpecialEKeys.AGE, 'prenatal')
                    return
                if unit == 'months':
                    number_part = number_part / 12
                    unit = 'years'

                if unit != 'years':
                    raise ValueError(f'Unexpected unit {unit} in age fuzzier')

                if number_part >= 80:
                    patch_meta.patch_value(SpecialEKeys.AGE, '80+')
                else:
                    base = math.floor(number_part / 10)
                    if base == 0:
                        patch_meta.patch_value(SpecialEKeys.AGE, '0-9')
                    else:
                        patch_meta.patch_value(SpecialEKeys.AGE, f'{base}0-{base}9')
