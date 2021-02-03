""" @see https://www.genenames.org/help/rest/

    Example JSON record at bottom of file
"""
import logging
from typing import Set, List

import requests

from annotation.models import CachedWebResource
from genes.gene_matching import GeneMatcher
from genes.models import HGNCImport, GeneAnnotationRelease, GeneSymbol, HGNC, GeneSymbolAlias
from genes.models_enums import HGNCStatus, GeneSymbolAliasSource
from library.django_utils import get_model_fields, get_field_counts
from library.utils import invert_dict

HGNC_BASE_URL = "http://rest.genenames.org/fetch/status/"


def store_hgnc_from_web(cached_web_resource: CachedWebResource):
    headers = {'Accept': 'application/json'}

    existing_hgnc_ids = set(HGNC.objects.all().values_list("pk", flat=True))

    for hgnc_status in HGNCStatus:
        url = HGNC_BASE_URL + hgnc_status.label
        logging.info("Fetching HGNC of status: %s", hgnc_status.label)
        r = requests.get(url, headers=headers)
        data = r.json()
        records = data["response"]["docs"]
        save_hgnc_records(existing_hgnc_ids, records)

    # Make sure gene symbols are matched to genes in each release
    for release in GeneAnnotationRelease.objects.all():
        gm = GeneMatcher(release)
        gm.match_unmatched_in_hgnc_and_gene_lists()

    status_counts = get_field_counts(HGNC.objects.all(), "status")
    cached_web_resource.description = ", ".join([f"{hs.label}: {status_counts[hs]}" for hs in HGNCStatus])
    cached_web_resource.save()


def save_hgnc_records(existing_hgnc_ids: Set, records: List):
    hgnc_status_lookup = invert_dict(dict(HGNCStatus.choices))
    hgnc_gene_names_new = []
    hgnc_gene_names_update = []
    gene_symbols = []
    gene_symbol_aliases = []

    hgnc_import = HGNCImport.objects.create()

    def _join_list(lst):
        return ",".join([str(s) for s in lst])

    for record in records:
        hgnc_id = record['hgnc_id']
        hgnc_id = int(hgnc_id.replace("HGNC:", ""))

        gene_symbol_id = record['symbol'].upper()
        gene_symbols.append(GeneSymbol(symbol=gene_symbol_id))

        previous_symbols = record.get('prev_symbol', [])
        alias_symbols = record.get('alias_symbol', [])
        status = hgnc_status_lookup[record['status']]

        def _get_list(key):
            return _join_list(record.get(key, []))

        hgnc = HGNC(pk=hgnc_id,
                    alias_symbols=_join_list(alias_symbols),
                    approved_name=record['name'],
                    ccds_ids=_get_list('ccds_id'),
                    ensembl_gene_id=record.get("ensembl_gene_id"),
                    gene_group_ids=_get_list('gene_group_id'),
                    gene_groups=_get_list('gene_group'),
                    gene_symbol_id=gene_symbol_id,
                    hgnc_import=hgnc_import,
                    location=record.get("location"),
                    mgd_ids=_get_list('mgd_id'),
                    omim_ids=_get_list('omim_id'),
                    previous_symbols=_join_list(previous_symbols),
                    refseq_ids=_get_list('refseq_accession'),
                    rgd_ids=_get_list('rgd_id'),
                    status=status,
                    ucsc_ids=record.get('ucsc_id'))

        if hgnc_id in existing_hgnc_ids:
            hgnc_gene_names_update.append(hgnc)
        else:
            existing_hgnc_ids.add(hgnc_id)
            hgnc_gene_names_new.append(hgnc)

        for alias_list in [previous_symbols, alias_symbols]:
            for alias in alias_list:
                alias = alias.strip().upper()
                gene_symbol_aliases.append(GeneSymbolAlias(alias=alias,
                                                           gene_symbol_id=gene_symbol_id,
                                                           source=GeneSymbolAliasSource.HGNC))

    if gene_symbols:
        logging.info("Upserting %d gene symbols", len(gene_symbols))
        GeneSymbol.objects.bulk_create(gene_symbols, ignore_conflicts=True)

    if gene_symbol_aliases:
        logging.info("Upserting %d gene symbol aliases", len(gene_symbol_aliases))
        GeneSymbolAlias.objects.bulk_create(gene_symbol_aliases, ignore_conflicts=True)

    if hgnc_gene_names_new:
        logging.info("Creating %d new hgnc_gene_names", len(hgnc_gene_names_new))
        HGNC.objects.bulk_create(hgnc_gene_names_new)

    if hgnc_gene_names_update:
        logging.info("Updating %d hgnc_gene_names", len(hgnc_gene_names_update))
        fields = get_model_fields(HGNC, ignore_fields=["id"])
        HGNC.objects.bulk_update(hgnc_gene_names_update, fields, batch_size=2000)


# http://rest.genenames.org/fetch/status/Approved
# http://rest.genenames.org/fetch/status/%22Entry%20Withdrawn%22


"""
    {'hgnc_id': 'HGNC:4171',
     'symbol': 'GATA2',
     'name': 'GATA binding protein 2',
     'status': 'Approved',
     'locus_type': 'gene with protein product',
     'prev_name': ['GATA-binding protein 2'],
     'alias_symbol': ['NFE1B'],
     'location': '3q21.3',
     'date_approved_reserved': '1992-11-03T00:00:00Z',
     'date_modified': '2020-12-06T00:00:00Z',
     'date_name_changed': '2001-11-28T00:00:00Z',
     'ena': ['AF169253'],
     'entrez_id': '2624',
     'mgd_id': ['MGI:95662'],
     'cosmic': 'GATA2',
     'orphanet': 274222,
     'pubmed_id': [1714909],
     'refseq_accession': ['NM_032638'],
     'gene_group': ['GATA zinc finger domain containing'],
     'vega_id': 'OTTHUMG00000159689',
     'lsdb': ['LRG_295|http://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_295.xml'],
     'ensembl_gene_id': 'ENSG00000179348',
     'ccds_id': ['CCDS3049', 'CCDS46903'],
     'locus_group': 'protein-coding gene',
     'omim_id': ['137295'],
     'uniprot_ids': ['P23769'],
     'ucsc_id': 'uc003eko.3',
     'rgd_id': ['RGD:2664'],
     'gene_group_id': [82],
     'location_sortable': '03q21.3',
     'agr': 'HGNC:4171',
     'symbol_report_tag': ['Stable symbol'],
     'uuid': '0c467aa5-50e2-4f56-acef-ce46ab6dba11',
     '_version_': 1690196606957125633}
"""