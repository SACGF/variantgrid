from django import db
from django.db import connection
import csv
import logging
import os

from library.utils import single_quote

LOCI_HEADER = ['contig_id', 'position', 'ref_id']
VARIANTS_HEADER = ['locus_id', 'alt_id']
COHORT_GENOTYPE_HEADER = ['collection_id', 'variant_id', "filters",
                          'ref_count', 'het_count', 'hom_count', 'unk_count', 'samples_zygosity',
                          'samples_allele_depth', 'samples_allele_frequency', 'samples_read_depth',
                          'samples_genotype_quality', 'samples_phred_likelihood', 'samples_filters']
MODIFIED_IMPORTED_VARIANT_HEADER = ['import_info_id', 'variant_id',
                                    'old_multiallelic', 'old_variant', 'old_variant_formatted']

GENE_COVERAGE_HEADER = [
    "gene_coverage_collection_id", "gene_symbol_id", "transcript_id", "transcript_version_id",
    "original_gene_symbol", "original_transcript", "size", "min", "max", "mean", "std_dev",
    "percent_1x", "percent_2x", "percent_5x", "percent_10x", "percent_15x", "percent_20x", "percent_25x",
    "percent_30x", "percent_40x", "percent_50x", "percent_60x", "percent_80x", "percent_100x",
    "percent_150x", "percent_200x", "percent_250x",
]

GENE_COVERAGE_CANONICAL_TRANSCRIPT_HEADER = ["canonical_transcript_collection_id"] + GENE_COVERAGE_HEADER



def get_database_settings():
    """ returns host, database_name, user, password """

    db_settings = db.connections.databases['default']
    host = db_settings['HOST']
    database_name = db_settings['NAME']
    user = db_settings['USER']
    password = db_settings['PASSWORD']
    return host, database_name, user, password


def sql_copy_csv(input_filename, table_name, columns, delimiter=',', quote=None):
    logging.info("sql_copy_csv %s", input_filename)
    with open(input_filename, 'rb') as f:
        return sql_copy_csv_file(f, table_name, columns, delimiter, quote=quote)


def sql_copy_csv_file(f, table_name, columns, delimiter=',', quote=None):
    cursor = connection.cursor()
    try:
        if quote:
            if delimiter != ',':
                msg = f"Don't know how to do this (sql_copy_csv_file sep='{delimiter}', quote='{quote}'"
                raise ValueError(msg)

            columns = ','.join(columns)
            quote = single_quote(quote)
            sql = f"COPY {table_name} ({columns}) FROM STDIN CSV QUOTE {quote};"
            value = cursor.copy_expert(sql, f)
        else:
            value = cursor.copy_from(f,
                                     table_name,
                                     delimiter,
                                     null='',
                                     columns=columns)

        logging.debug("copy returned %s - affected %d rows", value, cursor.rowcount)

        affected = cursor.rowcount
        return affected
    except Exception as e:
        logging.error(e)
        logging.info("columns = '%s'", columns)
        raise e
    finally:
        cursor.close()


def loci_sql_copy_csv(input_filename):
    return sql_copy_csv(input_filename, "snpdb_locus", LOCI_HEADER)


def variants_sql_copy_csv(input_filename):
    return sql_copy_csv(input_filename, "snpdb_variant", VARIANTS_HEADER)


def cohort_genotype_sql_copy_csv(input_filename, table_name):
    return sql_copy_csv(input_filename, table_name, COHORT_GENOTYPE_HEADER)


def modified_imported_variant_sql_copy_csv(input_filename):
    return sql_copy_csv(input_filename, "upload_modifiedimportedvariant", MODIFIED_IMPORTED_VARIANT_HEADER, quote='"')


def gene_coverage_sql_copy_csv(input_filename, table_name):
    return sql_copy_csv(input_filename, table_name, GENE_COVERAGE_HEADER, quote='"')


def gene_coverage_canonical_transcript_sql_copy_csv(input_filename, table_name):
    return sql_copy_csv(input_filename, table_name, GENE_COVERAGE_CANONICAL_TRANSCRIPT_HEADER, quote='"')


def write_sql_copy_csv(data, filename, **kwargs):
    if os.path.exists(filename):
        msg = f"We don't want to overwrite '{filename}'"
        raise ValueError(msg)

    with open(filename, "w") as f:
        writer = csv.writer(f, **kwargs)
        for row in data:
            try:
                writer.writerow(row)
            except Exception as e:
                logging.error("CSV couldn't write '%s'", row)
                raise e
#            f.write(','.join([str(x) if x is not None else '' for x in row]) + '\n')

    if not os.path.exists(filename):
        msg = f"write_sql_copy_csv didn't write file '{filename}' what happened?"
        raise ValueError(msg)
