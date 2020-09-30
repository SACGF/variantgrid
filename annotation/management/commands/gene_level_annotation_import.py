#!/usr/bin/env python3

"""

Imports Frank's Ensembl gene level annotation data

"""

from django.core.management.base import BaseCommand

from annotation.models import EnsemblGeneAnnotationVersion
from annotation.models.models_enums import TranscriptStatus, GenomicStrand
from genes.models import Gene
from genes.models_enums import AnnotationConsortium
from library.django_utils.django_file_utils import get_import_processing_filename
from library.file_utils import file_md5sum
from library.pandas_utils import df_nan_to_none
from snpdb.models.models_genome import GenomeBuild
from upload.vcf.sql_copy_files import write_sql_copy_csv, sql_copy_csv
import os
import numpy as np
import pandas as pd


ENSEMBL_GENE_ANNOTATION_HEADER = ['gene_id',
                                  'hgnc_symbol',
                                  'external_gene_name',
                                  'hgnc_symbol_lower',
                                  'hgnc_name',
                                  'synonyms',
                                  'previous_symbols',
                                  'hgnc_chromosome',
                                  'gene_family_tag',
                                  'gene_family_description',
                                  'hgnc_id',
                                  'entrez_gene_id',
                                  'uniprot_id',
                                  'ucsc_id',
                                  'omim_id',
                                  'enzyme_ids',
                                  'ccds_ids',
                                  'rgd_id',
                                  'mgi_id',
                                  'rvis_percentile',
                                  'refseq_gene_summary',
                                  'function_from_uniprotkb',
                                  'pathway_from_uniprotkb',
                                  'tissue_specificity_from_uniprotkb',
                                  'phenotypes_from_ensembl',
                                  'omim_phenotypes',
                                  'gene_biotype',
                                  'status',
                                  'chromosome_name',
                                  'start_position',
                                  'end_position',
                                  'band',
                                  'strand',
                                  'percentage_gc_content',
                                  'transcript_count',
                                  'in_cancer_gene_census']

# These needs to be ints, and needs to allow None
ENSEMBL_GENE_ANNOTATION_DTYPE = {"entrez_gene_id": np.object_,
                                 "start_position": np.object_,
                                 "end_position": np.object_,
                                 "omim_id": np.object_}


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('--genome-build', required=True, help='GenomeBuild name')
        parser.add_argument('--ensembl-version', required=True, type=int)
        parser.add_argument('--replace', action='store_true', help='Replace existing')
        parser.add_argument('ensembl_gene_annotation_tsv')

    def handle(self, *args, **options):
        filename = options["ensembl_gene_annotation_tsv"]
        ensembl_version = options["ensembl_version"]
        build_name = options["genome_build"]
        replace = options["replace"]

        self.stdout.write("Hashing file...")
        md5_hash = file_md5sum(filename)
        genome_build = GenomeBuild.get_name_or_alias(build_name)
        defaults = {"ensembl_version": ensembl_version,
                    "genome_build": genome_build}
        annotation_version, created = EnsemblGeneAnnotationVersion.objects.get_or_create(md5_hash=md5_hash,
                                                                                         defaults=defaults)

        genes_qs = Gene.objects.filter(annotation_consortium=AnnotationConsortium.ENSEMBL,
                                       geneversion__genome_build=genome_build).distinct()
        ensembl_genes_for_build = set(genes_qs.values_list("identifier", flat=True))
        if not ensembl_genes_for_build:
            raise ValueError(f"No genes for build {genome_build.name} - you need to import them first!")

        if created or replace:
            annotation_version.filename = filename
            annotation_version.save()
        else:
            self.stderr.write(f"EnsemblGeneAnnotationVersion with hd5_hash='{md5_hash}' already exists as pk={annotation_version.pk}. Use --replace to force")
            exit(1)

        version_id = annotation_version.pk
        self.stdout.write("Reading gene file...")
        df = pd.read_csv(filename, sep='\t', index_col=None, dtype=ENSEMBL_GENE_ANNOTATION_DTYPE)
        df = df_nan_to_none(df)
        records = []

        transcript_status_lookup = {v: k for (k, v) in TranscriptStatus.CHOICES}

        self.stdout.write("Parsing rows...")
        unmatched_genes = 0
        for _, row in df.iterrows():
            ensembl_gene_id = row['ensembl_gene_id']
            if ensembl_gene_id not in ensembl_genes_for_build:
                self.stdout.write(f"Couldn't find ensembl ID {ensembl_gene_id} in {genome_build}")
                unmatched_genes += 1
                continue

            row["gene_id"] = ensembl_gene_id
            row_data = row.reindex(ENSEMBL_GENE_ANNOTATION_HEADER)

            # issue #861 - Type conversions from Frank's spreadsheet to our new format
            row_data["status"] = transcript_status_lookup.get(row_data["status"])
            strand = row_data["strand"]  # Strand -1/1
            if float(strand) > 0:
                strand = GenomicStrand.SENSE
            else:
                strand = GenomicStrand.ANTISENSE
            row_data["strand"] = strand

            data = (version_id,) + tuple(row_data)
            records.append(data)

        self.stdout.write(f"Total: {len(records)} - unmatched: {unmatched_genes}")

        separator = '\t'  # Was having trouble with quoting CSVs
        self.stdout.write("Creating CSV to insert into database...\n")
        csv_filename = get_import_processing_filename(version_id, "ensembl_gene_annotation.csv", prefix='ensembl_gene_annotation')
        if os.path.exists(csv_filename):
            os.remove(csv_filename)
        write_sql_copy_csv(records, csv_filename, separator=separator)
        partition_table = annotation_version.get_partition_table()
        self.stdout.write(f"Inserting file '{csv_filename}' into partition {partition_table}\n")
        sql_copy_csv(csv_filename, partition_table, ['version_id'] + ENSEMBL_GENE_ANNOTATION_HEADER, separator=separator)
        self.stdout.write("Done!\n")
