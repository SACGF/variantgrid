"""
Imports HumanProteinAtlas data @see http://v15.proteinatlas.org/about/download
"""

import logging
import os

import pandas as pd
from django.core.management.base import BaseCommand, CommandError

from annotation.models import HumanProteinAtlasAnnotationVersion, HumanProteinAtlasTissueSample
from genes.models import GeneSymbol, Gene
from genes.models_enums import AnnotationConsortium
from library.django_utils.django_file_utils import get_import_processing_filename
from library.file_utils import file_md5sum
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv


def get_or_create_hpa_samples_ids(df) -> dict:
    """ returns a dict of name: pk """
    hpa_sample_ids = dict(HumanProteinAtlasTissueSample.objects.all().values_list("name", "pk"))
    new_records = []
    for hpa_sample_name in df['Tissue'].unique():
        if hpa_sample_name not in hpa_sample_ids:
            new_records.append(HumanProteinAtlasTissueSample(name=hpa_sample_name))
    if new_records:
        hpa_ts_records = HumanProteinAtlasTissueSample.objects.bulk_create(new_records, batch_size=2000)
        hpa_sample_ids.update({hpa_ts.name: hpa_ts.pk for hpa_ts in hpa_ts_records})
    return hpa_sample_ids


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--replace', action='store_true')
        parser.add_argument('--hpa-version', type=float, required=True)
        parser.add_argument('rna_tissue_tsv_zip')

    def handle(self, *args, **options):
        filename = options["rna_tissue_tsv_zip"]
        if not filename.endswith(".tsv.zip"):
            raise CommandError("Unknown file type, expecting 'rna_tissue.tsv.zip'")

        hpa_version_number = options["hpa_version"]
        if hpa_version_number < 21:
            raise CommandError("Only support v21 and later (untested as yet...)")

        logging.info("Hashing file...")
        md5_hash = file_md5sum(filename)
        hpa_version, created = HumanProteinAtlasAnnotationVersion.objects.get_or_create(md5_hash=md5_hash,
                                                                                        hpa_version=hpa_version_number,
                                                                                        defaults={"unit": "nTPM"})
        if created:
            hpa_version.filename = filename
            hpa_version.save()
        else:
            if options["replace"]:
                logging.warning("--replace option used, deleting and recreating partitions")
                hpa_version.delete_related_objects()
                hpa_version.create_partition()
            else:
                message = f"HumanProteinAtlasAnnotationVersion with hd5_hash='{md5_hash}' already exists as " \
                          f"pk={hpa_version.pk}. If you really want to import this again, use the --replace option"
                raise ValueError(message)

        version_id = hpa_version.pk
        logging.info("Reading hpa_version file...")
        df = pd.read_csv(filename, sep='\t', index_col=None)
        records = []

        hpa_samples_ids = get_or_create_hpa_samples_ids(df)
        ensembl_genes_qs = Gene.objects.filter(annotation_consortium=AnnotationConsortium.ENSEMBL)
        ensembl_genes = set(ensembl_genes_qs.values_list("pk", flat=True))
        gene_symbols = GeneSymbol.get_upper_case_lookup()
        new_gene_symbols = set()

        logging.info("Parsing rows...")
        for _, row in df.iterrows():
            gene_symbol_str = row['Gene name']
            ensembl_gene_str = row['Gene']
            sample_name = row['Tissue']
            value = row['nTPM']

            gene_symbol_id = gene_symbols.get(gene_symbol_str.upper())
            if gene_symbol_id is None:
                new_gene_symbols.add(GeneSymbol(symbol=gene_symbol_str))
                gene_symbol_id = gene_symbol_str

            if ensembl_gene_str in ensembl_genes:
                ensembl_gene_id = ensembl_gene_str
            else:
                ensembl_gene_id = None

            tissue_sample = hpa_samples_ids[sample_name]
            data = (version_id, gene_symbol_id, ensembl_gene_id, tissue_sample, value)
            records.append(data)

        if new_gene_symbols:
            GeneSymbol.objects.bulk_create(new_gene_symbols, batch_size=2000)

        delimiter = '\t'  # Was having trouble with quoting CSVs
        logging.info("Creating TSV to insert into database...")
        csv_filename = get_import_processing_filename(version_id, "human_protein_annotation.csv",
                                                      prefix='human_protein_atlas')
        if os.path.exists(csv_filename):
            os.remove(csv_filename)
        write_sql_copy_csv(records, csv_filename, delimiter=delimiter)
        partition_table = hpa_version.get_partition_table()
        logging.info("Inserting file '%s' into partition %s", csv_filename, partition_table)

        HPA_HEADER = ['version_id',
                      'gene_symbol_id',
                      'gene_id',
                      'tissue_sample_id',
                      'value']
        sql_copy_csv(csv_filename, partition_table, HPA_HEADER, delimiter=delimiter)
        logging.info("Done!")
