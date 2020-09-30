from django.db.models.aggregates import Count
from guardian.shortcuts import assign_perm
import logging
from toposort import toposort

import pandas as pd
from pedigree.models import Sex, PedFile, PedFileRecord, PedFileFamily, \
    create_automatch_pedigree
from snpdb.models import ImportStatus, Cohort


# @see http://www.sph.umich.edu/csg/abecasis/Pedstats/tour/input.html
# From https://genome.sph.umich.edu/wiki/Vt#Pedigree_File:

# FIELD            DESCRIPTION                                  VALID VALUES                           MISSING_VALUES
# Family ID        ID of this family                            [A-Za-z0-9_]+                          0
# Individual ID    ID(s) of this individual (comma separated)   [A-Za-z0-9_]+(,[A-Za-z0-9_]+)*         cannot be missing
# Paternal ID      ID of the father                             [A-Za-z0-9_]+                          0
# Maternal ID      ID of the mother                             [A-Za-z0-9_]+                          0
# Sex              Sex of the individual                        1=male, 2=female, other, male, female  other
# Phenotype        Phenotype                                    [A-Za-z0-9_]+                          -9

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


def save_ped_records(ped_file_family, family_df, dependency_graph):
    ped_records_dict = {}
    for samples in toposort(dependency_graph):
        for sample in samples:
            record = family_df.loc[sample]
            father_id = get_parent_id(record['father'])
            if father_id:
                father = ped_records_dict[father_id]
            else:
                father = None
            mother_id = get_parent_id(record['mother'])
            if mother_id:
                mother = ped_records_dict[mother_id]
            else:
                mother = None
            sex = get_sex(record['sex'])
            affection = get_affection(record['affection'])
            ped_record = PedFileRecord(family=ped_file_family,
                                       sample=sample,
                                       father=father,
                                       mother=mother,
                                       sex=sex,
                                       affection=affection)
            ped_record.save()
            ped_records_dict[sample] = ped_record


def create_ped_file_family(ped_file, family, family_df):
    dependency_graph = {}
    for sample_id, data in family_df.iterrows():
        dependency_graph[sample_id] = set()
        for parent in ['father', 'mother']:
            parent_id = get_parent_id(data[parent])
            if parent_id:
                dependency_graph[sample_id].add(parent_id)

    ped_family = PedFileFamily(ped_file=ped_file, name=family)
    ped_family.save()

    save_ped_records(ped_family, family_df, dependency_graph)
    return ped_family


def import_ped(ped_file, name, user):
    df = pd.read_csv(ped_file, sep=r'\s+', header=None, index_col=[0, 1],
                     usecols=range(6), names=PED_COLUMNS)
    ped_file = PedFile(user=user, name=name)
    ped_file.save()

    perm = 'pedigree.view_pedfile'
    assign_perm(perm, user, ped_file)
    for group in user.groups.all():
        assign_perm(perm, group, ped_file)

    families = []
    for family in df.index.levels[0]:
        family_df = df.loc[family]

        try:
            family = create_ped_file_family(ped_file, family, family_df)
            families.append(family)
            if family.errors:
                raise ValueError(family.errors)
            ped_file.import_status = ImportStatus.SUCCESS
        except Exception as e:
            logging.error("Exception: %s", e)
            ped_file.import_status = ImportStatus.ERROR
            raise
        finally:
            ped_file.save()

    return ped_file, families


def automatch_pedigree_samples(user, families, min_matching_samples):
    """ Create a pedigree if any sample names match user's cohort samples  """

    cohort_qs = Cohort.filter_for_user(user)

    for ped_file_family in families:
        sample_names = ped_file_family.pedfilerecord_set.values_list("sample", flat=True)
        qs = cohort_qs.filter(cohortsample__sample__name__in=sample_names)
        qs = qs.annotate(count=Count("id")).filter(count__gte=min_matching_samples)
        for cohort in qs:
            create_automatch_pedigree(user, ped_file_family, cohort)
