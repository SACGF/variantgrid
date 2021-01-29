#!/usr/bin/env python3

from django.core.management.base import BaseCommand
import logging

from genes.gene_matching import GeneMatcher
from genes.models import HGNCStatus, HGNCGeneNamesImport, HGNC, GeneSymbolAlias, GeneSymbol, \
    GeneAnnotationRelease
from genes.models_enums import GeneSymbolAliasSource
from library.django_utils import get_model_fields
from library.pandas_utils import df_replace_nan
from library.utils import invert_dict
import pandas as pd

HGNC_ID = "HGNC ID"
APPROVED_SYMBOL = "Approved symbol"
APPROVED_NAME = "Approved name"
STATUS = "Status"
PREVIOUS_SYMBOLS = "Previous symbols"
ALIAS_SYMBOLS = "Alias symbols"
REFSEQ_IDS = "RefSeq IDs"


def insert_hgnc_data(hgnc_import, existing_hgnc_ids: set, df):
    hgnc_status_lookup = invert_dict(dict(HGNCStatus.choices))
    hgnc_gene_names_new = []
    hgnc_gene_names_update = []
    gene_symbols = []
    gene_symbol_aliases = []

    for _, row in df.iterrows():
        hgnc_id = row[HGNC_ID]
        hgnc_id = int(hgnc_id.replace("HGNC:", ""))

        gene_symbol_id = row[APPROVED_SYMBOL].upper()
        gene_symbols.append(GeneSymbol(symbol=gene_symbol_id))

        status = hgnc_status_lookup[row[STATUS]]
        hgnc = HGNC(pk=hgnc_id,
                    hgnc_import=hgnc_import,
                    gene_symbol_id=gene_symbol_id,
                    approved_name=row[APPROVED_NAME],
                    status=status,
                    previous_symbols=row[PREVIOUS_SYMBOLS],
                    alias_symbols=row[ALIAS_SYMBOLS],
                    refseq_ids=row[REFSEQ_IDS])
        if hgnc_id in existing_hgnc_ids:
            hgnc_gene_names_update.append(hgnc)
        else:
            hgnc_gene_names_new.append(hgnc)

        for gene_column in [PREVIOUS_SYMBOLS, ALIAS_SYMBOLS]:
            genes_str = row[gene_column]
            if genes_str:
                for alias in genes_str.split(","):
                    alias = alias.strip().upper()
                    gene_symbol_aliases.append(GeneSymbolAlias(alias=alias,
                                                               gene_symbol_id=gene_symbol_id,
                                                               source=GeneSymbolAliasSource.HGNC))

    if gene_symbols:
        GeneSymbol.objects.bulk_create(gene_symbols, ignore_conflicts=True)

    if gene_symbol_aliases:
        logging.info("Creating %d gene symbol aliases", len(gene_symbol_aliases))
        GeneSymbolAlias.objects.bulk_create(gene_symbol_aliases, ignore_conflicts=True)

    if hgnc_gene_names_new:
        logging.info("Creating %d new hgnc_gene_names", len(hgnc_gene_names_new))
        HGNC.objects.bulk_create(hgnc_gene_names_new)

    if hgnc_gene_names_update:
        logging.info("Updating %d hgnc_gene_names", len(hgnc_gene_names_update))
        fields = get_model_fields(HGNC, ignore_fields=["id"])
        HGNC.objects.bulk_update(hgnc_gene_names_update, fields, batch_size=2000)


class Command(BaseCommand):

    def add_arguments(self, parser):
        parser.add_argument('hgnc_gene_symbols_tsv')

    def handle(self, *args, **options):
        hgnc_gene_symbols_tsv = options["hgnc_gene_symbols_tsv"]

        hgnc_import = HGNCGeneNamesImport.objects.create()

        df = pd.read_csv(hgnc_gene_symbols_tsv, sep='\t', index_col=False)
        df = df_replace_nan(df)

        for c in [HGNC_ID, APPROVED_SYMBOL, APPROVED_NAME, STATUS, PREVIOUS_SYMBOLS, ALIAS_SYMBOLS, REFSEQ_IDS]:
            if c not in df.columns:
                url = "https://www.genenames.org/cgi-bin/download"
                msg = f"TSV ({df.columns}) is missing required column '{c}', is it downloaded from '{url}'?"
                raise ValueError(msg)

        existing_hgnc_ids = set(HGNC.objects.all().values_list("pk", flat=True))
        insert_hgnc_data(hgnc_import, existing_hgnc_ids, df)

        # Make sure gene symbols are matched to genes in each release
        for release in GeneAnnotationRelease.objects.all():
            gm = GeneMatcher(release)
            gm.match_unmatched_in_hgnc_and_gene_lists()

