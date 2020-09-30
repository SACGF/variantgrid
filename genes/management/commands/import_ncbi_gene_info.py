#!/usr/bin/env python3

# Data can be downloaded here:
# * ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
import re

from django.core.management.base import BaseCommand
from django.db.models.functions import Upper

from genes.models import GeneSymbol, GeneSymbolAlias
import pandas as pd

from genes.models_enums import GeneSymbolAliasSource


class Command(BaseCommand):
    BULK_INSERT_SIZE = 2000

    def add_arguments(self, parser):
        parser.add_argument('gene_info', help="eg Homo_sapiens.gene_info.gz")

    def handle(self, *args, **options):
        gene_info_filename = options["gene_info"]
        known_gene_symbols = GeneSymbol.get_upper_case_lookup()
        known_symbols = set(known_gene_symbols) | set(GeneSymbolAlias.get_upper_case_lookup())

        GeneSymbolAlias.objects.filter(source=GeneSymbolAliasSource.NCBI).delete()
        gene_info_df = pd.read_csv(gene_info_filename, sep='\t')
        symbols_synonyms = gene_info_df[["Symbol", "Synonyms"]]
        symbols_synonyms.loc[:, "Symbol"] = symbols_synonyms.loc[:, "Symbol"].str.upper()
        symbols_synonyms.loc[:, "Synonyms"] = symbols_synonyms.loc[:, "Synonyms"].str.upper()

        gene_symbols = []
        gene_symbol_aliases = []
        for _, (symbol, synonyms) in symbols_synonyms.iterrows():
            gene_symbol_id = known_gene_symbols.get(symbol)
            if gene_symbol_id is None:
                gene_symbols.append(GeneSymbol(symbol=symbol))
                gene_symbol_id = symbol

            if synonyms == "-":
                continue

            for s in synonyms.split("|"):
                if s not in known_symbols:
                    known_symbols.add(s)
                    gene_symbol_aliases.append(GeneSymbolAlias(alias=s,
                                                               gene_symbol_id=gene_symbol_id,
                                                               source=GeneSymbolAliasSource.NCBI))

        if gene_symbols:
            GeneSymbol.objects.bulk_create(gene_symbols, ignore_conflicts=True)
        if gene_symbol_aliases:
            print(f"Creating {len(gene_symbol_aliases)} gene symbol aliases")
            GeneSymbolAlias.objects.bulk_create(gene_symbol_aliases, ignore_conflicts=True)
