""" Retire GeneInfo (#693).

    GeneInfo tagged gene symbols with small icons on GeneGrid via 3 rows created in
    genes/migrations/0002_initial_data.py - Alternative Haplotype, Pseudogene and Triplet Repeat
    Disorders - each needing a gene list uploaded then hand-wired in admin before anything showed.

    Triplet repeats are covered by PanelApp STR entities (already on the gene page), alternative
    haplotypes were a GRCh37 artefact, and pseudogenes will be revisited via GIAB challenging
    medically relevant genes as a CachedWebResource.

    The "GeneInfo" GeneListCategory and any gene lists in it are left alone - on deployments that did
    the manual wiring those are real uploaded gene lists.
"""
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("genes", "0084_panel_app_local_cache_gene_symbol_hgnc"),
    ]

    operations = [
        migrations.DeleteModel(
            name="GeneInfo",
        ),
    ]
