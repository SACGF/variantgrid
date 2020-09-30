"""
Imports HumanProteinAtlas data @see http://v15.proteinatlas.org/about/download
"""

from django.core.management.base import BaseCommand, CommandError
import logging
import os

from annotation.models import HumanProteinAtlasAnnotationVersion, HumanProteinAtlasTissueSample, \
    HumanProteinAtlasAbundance
from library.django_utils.django_file_utils import get_import_processing_filename
from library.file_utils import file_md5sum
import pandas as pd
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv


def get_or_create_hpa_samples_ids(df):
    """ returns a dict of name: pk """
    hpa_sample_ids = {}

    for hpa_sample_name in df['Sample'].unique():
        hpa_sample, _ = HumanProteinAtlasTissueSample.objects.get_or_create(name=hpa_sample_name)
        hpa_sample_ids[hpa_sample_name] = hpa_sample.pk

    return hpa_sample_ids


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--replace', action='store_true')
        parser.add_argument('--hpa-version', type=float, required=True)
        parser.add_argument('rna_tissue_csv_zip')

    def handle(self, *args, **options):
        filename = options["rna_tissue_csv_zip"]
        if not filename.endswith(".csv.zip"):
            raise CommandError("Unknown file type, expecting 'rna_tissue.csv.zip'")

        hpa_version_number = options["hpa_version"]
        logging.info("Hashing file...")
        md5_hash = file_md5sum(filename)
        hpa_version, created = HumanProteinAtlasAnnotationVersion.objects.get_or_create(md5_hash=md5_hash,
                                                                                        hpa_version=hpa_version_number)
        if created:
            hpa_version.filename = filename
            hpa_version.save()
        else:
            if options["replace"]:
                logging.warning("--replace option used, deleting and recreating partitions")
                hpa_version.delete_related_objects()
                hpa_version.create_partition()
            else:
                logging.error("HumanProteinAtlasAnnotationVersion with hd5_hash='%s' already exists as pk=%d.", md5_hash, hpa_version.pk)
                logging.error("If you really want to import this again, use the --replace option")
                exit(1)

        version_id = hpa_version.pk
        logging.info("Reading hpa_version file...")
        df = pd.read_csv(filename, index_col=None)
        records = []

        hpa_samples_ids = get_or_create_hpa_samples_ids(df)
        abundance_lookup = {b: a for (a, b) in HumanProteinAtlasAbundance.CHOICES}

        logging.info("Parsing rows...")
        for (_, row) in df.iterrows():
            ensembl_gene_id = row['Gene']
            sample_name = row['Sample']
            abundance_string = row['Abundance']

            tissue_sample = hpa_samples_ids[sample_name]
            abundance = abundance_lookup[abundance_string]
            data = (version_id, ensembl_gene_id, tissue_sample, abundance)
            records.append(data)

        separator = '\t'  # Was having trouble with quoting CSVs
        logging.info("Creating TSV to insert into database...")
        csv_filename = get_import_processing_filename(version_id, "human_protein_annotation.csv", prefix='ensembl_gene_annotation')
        if os.path.exists(csv_filename):
            os.remove(csv_filename)
        write_sql_copy_csv(records, csv_filename, separator=separator)
        partition_table = hpa_version.get_partition_table()
        logging.info("Inserting file '%s' into partition %s", csv_filename, partition_table)

        HPA_HEADER = ['version_id',
                      'gene_id',
                      'tissue_sample_id',
                      'abundance']
        sql_copy_csv(csv_filename, partition_table, HPA_HEADER, separator=separator)
        logging.info("Done!")
