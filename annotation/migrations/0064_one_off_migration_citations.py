# Generated by Django 4.0.7 on 2022-12-06 03:14

import json

from django.db import migrations


def _migrate_citations(apps, _schema_editor):
    print("About to populate Citation2")
    Citation = apps.get_model("annotation", "Citation")
    Citation2 = apps.get_model("annotation", "Citation2")
    CachedCitation = apps.get_model("annotation", "CachedCitation")
    all_citation_2s = []

    citation_map = {}

    for citation in Citation.objects.all():
        if len(all_citation_2s) % 1000 == 0:
            print(f"Migrated {len(all_citation_2s)} citations")

        citation_json = None
        if cached := CachedCitation.objects.filter(citation=citation).first():
            if json_string := cached.json_string:
                try:
                    citation_json = json.loads(json_string)
                except RuntimeError as re:
                    print(re)

        source = citation.citation_source
        dest_source = None
        dest_index = f"{citation.citation_id}"

        if source == "P":
            dest_source = "PMID"

        elif source == "C":
            dest_source = "PMCID"
            if not dest_index.startswith("PMC"):
                dest_index = f"PMC{dest_index}"

        elif source == "N":
            dest_source = "BookShelf"
            if not dest_index.startswith("NBK"):
                dest_index = f"NBK{dest_index}"

        dest_id = f"{dest_source}:{dest_index}"
        all_citation_2s.append(Citation2(
            pk=dest_id,
            old_id=citation.pk,
            source=dest_source,
            index=dest_index,
            data_json=citation_json
        ))
        citation_map[citation.pk] = dest_id

    Citation2.objects.bulk_create(objs=all_citation_2s, ignore_conflicts=True, batch_size=1000)
    print(f"Completed {len(all_citation_2s)} migrations")

    all_gene_symbol_citations = []
    GeneSymbolCitation = apps.get_model("annotation", "GeneSymbolCitation")
    for gene_symbol_citation in GeneSymbolCitation.objects.all():
        if len(all_gene_symbol_citations) % 1000 == 0:
            print(f"Migrated {len(all_gene_symbol_citations)} gene_symbol_citation")
        old_id = gene_symbol_citation.citation_id
        gene_symbol_citation.citation2_id = citation_map[old_id]
        all_gene_symbol_citations.append(gene_symbol_citation)
    print(f"About to bulk update {len(all_gene_symbol_citations)} gene_symbol_citations - this may take a few minutes")
    GeneSymbolCitation.objects.bulk_update(objs=all_gene_symbol_citations, fields=['citation2'], batch_size=1000)
    print(f"Completed {len(all_gene_symbol_citations)} migrations")

    all_clinvar_citations = []
    ClinVarCitation = apps.get_model("annotation", "ClinVarCitation")
    for clinvar_citation in ClinVarCitation.objects.filter(citation__isnull=False):
        old_id = clinvar_citation.citation_id
        clinvar_citation.citation2_id = citation_map[old_id]
        all_clinvar_citations.append(clinvar_citation)
    print(f"About to bulk update {len(all_clinvar_citations)} clinvar_citations - this may take a few minutes")
    ClinVarCitation.objects.bulk_update(objs=all_clinvar_citations, fields=['citation2'], batch_size=10000)
    print(f"Completed {len(all_clinvar_citations)} migrations")


def _unmigrate_citations(apps, _schema_editor):
    Citation2 = apps.get_model("annotation", "Citation2")
    GeneSymbolCitation = apps.get_model("annotation", "GeneSymbolCitation")
    ClinVarCitation = apps.get_model("annotation", "ClinVarCitation")

    ClinVarCitation.objects.all().update(citation2=None)
    GeneSymbolCitation.objects.all().update(citation2=None)
    Citation2.objects.all().delete()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0063_citation2_clinvarcitation_citation2_and_more'),
    ]

    operations = [
        migrations.RunPython(_migrate_citations, _unmigrate_citations)
    ]
