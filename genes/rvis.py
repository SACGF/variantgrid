from io import BytesIO

import pandas as pd
import requests

from annotation.models import CachedWebResource
from genes.models import RVIS, GeneSymbol


def store_rvis_from_web(cached_web_resource: CachedWebResource):
    RVIS_URL = "http://genic-intolerance.org/data/GenicIntolerance_v3_12Mar16.txt"

    r = requests.get(RVIS_URL)
    f = BytesIO(r.content)
    save_rvis_records(cached_web_resource, f)

    cached_web_resource.description = f"{RVIS.objects.count()} RVIS records"
    cached_web_resource.save()


def save_rvis_records(cached_web_resource: CachedWebResource, f):
    GENE_COLUMN = "GENE"
    OE_RATIO_PERCENTILE = "OEratio-percentile[ExAC]"

    known_gene_symbols = GeneSymbol.get_upper_case_lookup()

    df = pd.read_csv(f, sep='\t')
    gene_symbols = []
    rvis_records = []
    for _, row in df[[GENE_COLUMN, OE_RATIO_PERCENTILE]].iterrows():
        symbol = row[GENE_COLUMN]
        oe_ratio_percentile = row[OE_RATIO_PERCENTILE]

        # Grab existing (case insensitive)
        gene_symbol_id = known_gene_symbols.get(symbol.upper())
        if gene_symbol_id is None:
            # Will create new one
            gene_symbol_id = symbol
            gene_symbols.append(GeneSymbol(symbol=gene_symbol_id))

        rvis_records.append(RVIS(cached_web_resource=cached_web_resource, gene_symbol_id=gene_symbol_id,
                                 oe_ratio_percentile=oe_ratio_percentile))

    if gene_symbols:
        GeneSymbol.objects.bulk_create(gene_symbols)

    RVIS.objects.bulk_create(rvis_records)
