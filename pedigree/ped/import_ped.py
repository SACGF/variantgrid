import logging
from graphlib import CycleError

import pandas as pd
from django.db import transaction
from django.db.models.aggregates import Count
from guardian.shortcuts import assign_perm

from library.utils.collection_utils import toposort_groups
from pedigree.models import (
    PedFile,
    PedFileFamily,
    PedFileRecord,
    Pedigree,
    create_automatch_pedigree,
)
from pedigree.ped.ped_file_utils import PED_COLUMNS, get_affection, get_parent_id, get_sex
from snpdb.models import Cohort, ImportStatus


def save_ped_records(ped_file_family, family_df, dependency_graph):
    ped_records_dict = {}
    for samples in toposort_groups(dependency_graph):
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
    """ Every family is validated before any is kept: an invalid one (no affected member, a parent
        that isn't a record, a parent cycle) rolls the families back, marks the PedFile ERROR and
        raises with the errors of every bad family """
    df = pd.read_csv(ped_file, sep=r'\s+', header=None, index_col=[0, 1],
                     usecols=range(6), names=PED_COLUMNS)
    ped_file = PedFile(user=user, name=name)
    ped_file.save()

    perm = 'pedigree.view_pedfile'
    assign_perm(perm, user, ped_file)
    for group in user.groups.all():
        assign_perm(perm, group, ped_file)

    try:
        families = _create_ped_file_families(ped_file, df)
        ped_file.import_status = ImportStatus.SUCCESS
    except Exception as e:
        logging.error("Exception: %s", e)
        ped_file.import_status = ImportStatus.ERROR
        raise
    finally:
        ped_file.save()

    return ped_file, families


def _create_ped_file_families(ped_file, df):
    family_errors = {}
    families = []
    with transaction.atomic():
        for family_name in df.index.levels[0]:
            family_df = df.loc[family_name]
            try:
                family = create_ped_file_family(ped_file, family_name, family_df)
            except (KeyError, CycleError) as e:
                family_errors[family_name] = [f"{e.__class__.__name__}: {e}"]
                continue
            if family.errors:
                family_errors[family_name] = family.errors
            families.append(family)
        if family_errors:
            raise ValueError(family_errors)
    return families


def automatch_pedigree_samples(user, families, min_matching_samples):
    """ Create a pedigree if any sample names match user's cohort samples. Re-uploading a PED file
        makes a new PedFile and families, so a cohort that already has a pedigree over a family of
        the same name and sample names is skipped """

    cohort_qs = Cohort.filter_for_user(user)

    for ped_file_family in families:
        sample_names = set(ped_file_family.pedfilerecord_set.values_list("sample", flat=True))
        qs = cohort_qs.filter(cohortsample__sample__name__in=sample_names)
        qs = qs.annotate(count=Count("id")).filter(count__gte=min_matching_samples)
        for cohort in qs:
            if not _cohort_has_equivalent_pedigree(cohort, ped_file_family.name, sample_names):
                create_automatch_pedigree(user, ped_file_family, cohort)


def _cohort_has_equivalent_pedigree(cohort, family_name, sample_names) -> bool:
    for pedigree in Pedigree.objects.filter(cohort=cohort, ped_file_family__name=family_name):
        existing = set(pedigree.ped_file_family.pedfilerecord_set.values_list("sample", flat=True))
        if existing == sample_names:
            return True
    return False
