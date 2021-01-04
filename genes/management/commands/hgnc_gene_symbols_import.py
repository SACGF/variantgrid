#!/usr/bin/env python3

"""
HGNC Gene Symbols
"""

from django.core.management.base import BaseCommand
import logging

from genes.models import HGNCStatus, HGNCGeneNamesImport, HGNCGeneNames, GeneSymbolAlias, GeneSymbol
from genes.models_enums import GeneSymbolAliasSource
from library.pandas_utils import df_replace_nan
from library.utils import invert_dict
import pandas as pd

HGNC_ID = "HGNC ID"
APPROVED_SYMBOL = "Approved symbol"
APPROVED_NAME = "Approved name"
STATUS = "Status"
PREVIOUS_SYMBOLS = "Previous symbols"
SYNONYMS = "Synonyms"
REFSEQ_IDS = "RefSeq IDs"


def insert_hgnc_data(hgnc_import, df):
    hgnc_status_lookup = invert_dict(dict(HGNCStatus.choices))
    hgnc_gene_names = []
    gene_symbols = []
    gene_symbol_aliases = []

    for _, row in df.iterrows():
        hgnc_id = row[HGNC_ID]
        hgnc_id = int(hgnc_id.replace("HGNC:", ""))

        gene_symbol_id = row[APPROVED_SYMBOL].upper()
        gene_symbols.append(GeneSymbol(symbol=gene_symbol_id))

        status = hgnc_status_lookup[row[STATUS]]
        hgnc = HGNCGeneNames(pk=hgnc_id,
                             hgnc_import=hgnc_import,
                             approved_symbol=gene_symbol_id,
                             approved_name=row[APPROVED_NAME],
                             status=status,
                             previous_symbols=row[PREVIOUS_SYMBOLS],
                             synonyms=row[SYNONYMS],
                             refseq_ids=row[REFSEQ_IDS])
        hgnc_gene_names.append(hgnc)

        for gene_column in [PREVIOUS_SYMBOLS, SYNONYMS]:
            genes_str = row[gene_column]
            if genes_str:
                for alias in genes_str.split(","):
                    alias = alias.strip().upper()
                    gene_symbol_aliases.append(GeneSymbolAlias(alias=alias,
                                                               gene_symbol_id=gene_symbol_id,
                                                               source=GeneSymbolAliasSource.HGNC))

    if gene_symbols:
        GeneSymbol.objects.bulk_create(gene_symbols, ignore_conflicts=True)
    if hgnc_gene_names:
        logging.info("Creating %d hgnc_gene_names", len(hgnc_gene_names))
        HGNCGeneNames.objects.bulk_create(hgnc_gene_names)
    if gene_symbol_aliases:
        logging.info("Creating %d gene symbol aliases", len(gene_symbol_aliases))
        GeneSymbolAlias.objects.bulk_create(gene_symbol_aliases, ignore_conflicts=True)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('hgnc_gene_symbols_tsv')

    def handle(self, *args, **options):
        hgnc_gene_symbols_tsv = options["hgnc_gene_symbols_tsv"]

        HGNCGeneNamesImport.objects.all().delete()  # Only want 1
        hgnc_import = HGNCGeneNamesImport.objects.create()

        df = pd.read_csv(hgnc_gene_symbols_tsv, sep='\t', index_col=False)
        df = df_replace_nan(df)

        REQUIRED_COLUMNS = [HGNC_ID, APPROVED_SYMBOL, APPROVED_NAME, STATUS, PREVIOUS_SYMBOLS, SYNONYMS, REFSEQ_IDS]
        for c in REQUIRED_COLUMNS:
            if c not in df.columns:
                url = "https://www.genenames.org/cgi-bin/download"
                msg = f"TSV ({df.columns}) is missing required column '{c}', is it downloaded from '{url}'?"
                raise ValueError(msg)

        insert_hgnc_data(hgnc_import, df)
