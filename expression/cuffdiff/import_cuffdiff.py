import logging
import os

import pandas as pd
from django.conf import settings
from django.db import connection
from guardian.shortcuts import assign_perm

from expression.models import CuffDiffFile, CuffDiffRecord
from library.pandas_utils import df_handle_below_minimum_floats
from snpdb.models.models_enums import ImportStatus, AnnotationLevel

# None means don't validate
COLUMNS = ['test_id', None, 'gene', 'locus', None, None, 'status', 'value_1', 'value_2', 'log2(fold_change)', 'test_stat', 'p_value', 'q_value', 'significant']
BULK_INSERT_SIZE = 5000


def validate_columns(df, annotation_level):
    expected_cols = COLUMNS[1:]  # Index not counted
    if annotation_level == AnnotationLevel.GENE:
        expected_cols[0] = 'gene_id'
    elif annotation_level == AnnotationLevel.TRANSCRIPT:
        expected_cols[0] = 'transcript_id'
    else:
        raise ValueError(f"Bad value for annotation level {annotation_level}")

    num_expected_cols = len(expected_cols)
    num_actual_cols = len(df.columns)
    if num_expected_cols != num_actual_cols:
        msg = f"Expected Cuffdiff file ({annotation_level}) to have {num_expected_cols} columns, got {num_actual_cols}!"
        raise ValueError(msg)

    for i, (expected, actual) in enumerate(zip(expected_cols, df.columns)):
        if expected:
            if expected != actual:
                msg = f"Expected Cuffdiff file column {i} to be '{expected}', was: '{actual}'"
                raise ValueError(msg)


def insert_cuffdiff_records(cuff_diff_file, df):
    inserted = 0

    records = []
    for test_id, row in df.iterrows():
        cdr = CuffDiffRecord(cuff_diff_file=cuff_diff_file,
                             test_id=test_id,
                             reference_id=row[0],
                             locus=row['locus'],
                             status=row['status'],
                             value_1=row['value_1'],
                             value_2=row['value_2'],
                             log2_fold_change=row['log2(fold_change)'],
                             test_stat=row['test_stat'],
                             p_value=row['p_value'],
                             q_value=row['q_value'],)
        records.append(cdr)
        num_records = len(records)
        if num_records >= BULK_INSERT_SIZE:
            inserted += num_records
            logging.info("Inserting cuffdiff records")
            CuffDiffRecord.objects.bulk_create(records)
            records = []

    num_records = len(records)
    if num_records:
        inserted += num_records
        CuffDiffRecord.objects.bulk_create(records)

    return num_records


def get_samples_dict(df):
    """ Returns dict of {'sample_1' : name1, 'sample_2 : 'name2'} """
    samples = {}
    for sample_col in ['sample_1', 'sample_2']:
        sample_set = set(df[sample_col])
        if len(sample_set) != 1:
            msg = f"Can only handle cuffdiff files with 2 labels. '{sample_col}' should have a single value, has: '{sample_set}'"
            raise ValueError(msg)
        sample, = sample_set
        samples[sample_col] = sample
    return samples


def link_records(cuff_diff_file):
    EXPRESSION_SCRIPTS_DIR = os.path.join(settings.BASE_DIR, "expression", "dbscripts", settings.BACKEND_ENGINE)

    if cuff_diff_file.annotation_level == AnnotationLevel.GENE:
        script_name = 'link_cuffdiff_gene_ids.sql'
    elif cuff_diff_file.annotation_level == AnnotationLevel.TRANSCRIPT:
        script_name = 'link_cuffdiff_transcript_ids.sql'
    else:
        msg = f"Unknown annotation level '{cuff_diff_file.annotation_level}'"
        raise ValueError(msg)

    with open(os.path.join(EXPRESSION_SCRIPTS_DIR, script_name)) as f:
        raw_sql = f.read()
    sql = raw_sql % {"cuff_diff_file_id": cuff_diff_file.pk}

    cursor = connection.cursor()
    try:
        cursor.execute(sql)
        affected_rows = cursor.rowcount
        logging.info("affected %d rows", affected_rows)
        return affected_rows
    finally:
        cursor.close()


def get_annotation_level_from_columns(df):
    try:
        col_name = df.columns[0]  # Index doesn't count
        if col_name == 'gene_id':
            annotation_level = AnnotationLevel.GENE
        elif col_name == 'transcript_id':
            annotation_level = AnnotationLevel.TRANSCRIPT
        else:
            msg = f"Unknown column 1 value of '{col_name}'"
            raise Exception(msg)
    except Exception as e:
        msg = f"Could not auto-detect annotation level: {e}"
        raise ValueError(msg)

    return annotation_level


def import_cuffdiff(cuffdiff_file, name, user, annotation_level=None):
    logging.debug("import_cuffdiff(%s)", name)

    df = pd.read_csv(cuffdiff_file, sep='\t')
    logging.debug("columns: %s", df.columns)
    # Fix for Postgres "out of range for type double precision" error
    df = df_handle_below_minimum_floats(df)

    if annotation_level is None:
        annotation_level = get_annotation_level_from_columns(df)

    validate_columns(df, annotation_level)
    samples = get_samples_dict(df)
    logging.debug("samples = %s", samples)

    cuff_diff_file = CuffDiffFile(user=user,
                                  name=name,
                                  annotation_level=annotation_level,
                                  import_status=ImportStatus.IMPORTING,
                                  sample_1=samples['sample_1'],
                                  sample_2=samples['sample_2'],)
    cuff_diff_file.save()

    # Allow others in group to see it
    perm = 'expression.view_cuff_diff_file'
    assign_perm(perm, user, cuff_diff_file)
    for group in user.groups.all():
        assign_perm(perm, group, cuff_diff_file)

    imported_records = insert_cuffdiff_records(cuff_diff_file, df)

    cuff_diff_file.imported_records = imported_records
    cuff_diff_file.matched_records = link_records(cuff_diff_file)
    cuff_diff_file.import_status = ImportStatus.SUCCESS
    cuff_diff_file.save()

    return cuff_diff_file
