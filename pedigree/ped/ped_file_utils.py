# @see http://www.sph.umich.edu/csg/abecasis/Pedstats/tour/input.html
# From https://genome.sph.umich.edu/wiki/Vt#Pedigree_File:

# FIELD            DESCRIPTION                                  VALID VALUES                           MISSING_VALUES
# Family ID        ID of this family                            [A-Za-z0-9_]+                          0
# Individual ID    ID(s) of this individual (comma separated)   [A-Za-z0-9_]+(,[A-Za-z0-9_]+)*         cannot be missing
# Paternal ID      ID of the father                             [A-Za-z0-9_]+                          0
# Maternal ID      ID of the mother                             [A-Za-z0-9_]+                          0
# Sex              Sex of the individual                        1=male, 2=female, other, male, female  other
# Phenotype        Phenotype                                    [A-Za-z0-9_]+                          -9
from patients.models_enums import Sex

PED_COLUMNS = ['family', 'sample', 'father', 'mother', 'sex', 'affection']

SEX_LOOKUP = {'1': Sex.MALE,
              'M': Sex.MALE,
              '2': Sex.FEMALE,
              'F': Sex.FEMALE}

UNKNOWN_PARENT_VALUES = ['.', 0, 'Unknown', '0']
AFFECTION_VALUES = {'X': None,
                    '0': None,
                    'U': False,
                    '1': False,
                    'A': True,
                    '2': True,
                    '-9': None}  # Phenotips uses this for unknown values


def get_parent_id(parent_id):
    if parent_id in UNKNOWN_PARENT_VALUES:
        return None
    return parent_id


def get_sex(sex):
    sex = str(sex)
    return SEX_LOOKUP.get(sex, None)


def get_affection(affection):
    affection = str(affection)
    if affection not in AFFECTION_VALUES:
        msg = f"Affection '{affection}' unknown (one of: {AFFECTION_VALUES})"
        raise ValueError(msg)
    return AFFECTION_VALUES[affection]