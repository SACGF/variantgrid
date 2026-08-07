""" Key cached PanelApp genes off HGNC ID rather than gene symbol (#1667).

    PanelApp reports gene symbols from a dated Ensembl/HGNC snapshot, so the symbol we cached could
    lag - or mean a different gene than - our current approved symbol. Their records carry an HGNC ID,
    which is stable, so PanelAppPanelLocalCacheGeneSymbol now points at HGNC and keeps their symbol
    as gene_symbol_reported for provenance.

    Existing rows carry forward: gene_symbol_reported comes from the old gene_symbol FK (whose pk is
    the symbol string), and hgnc is resolved from the stored API response. Rows whose HGNC ID we have
    no record of, or that predate PanelApp returning one, keep hgnc=None and fall back to the
    reported symbol.
"""
from django.db import migrations, models
import django.db.models.deletion


def backfill_hgnc_and_reported_symbol(apps, schema_editor):
    PanelAppPanelLocalCacheGeneSymbol = apps.get_model("genes", "PanelAppPanelLocalCacheGeneSymbol")
    HGNC = apps.get_model("genes", "HGNC")

    known_hgnc_pks = set(HGNC.objects.all().values_list("pk", flat=True))
    to_update = []
    for record in PanelAppPanelLocalCacheGeneSymbol.objects.all().iterator(chunk_size=2000):
        record.gene_symbol_reported = record.gene_symbol_id or ""

        hgnc_pk = None
        hgnc_id = ((record.data or {}).get("gene_data") or {}).get("hgnc_id")
        if hgnc_id:
            try:
                hgnc_pk = int(str(hgnc_id).replace("HGNC:", ""))
            except ValueError:
                hgnc_pk = None
        record.hgnc_id = hgnc_pk if hgnc_pk in known_hgnc_pks else None

        to_update.append(record)
        if len(to_update) >= 2000:
            PanelAppPanelLocalCacheGeneSymbol.objects.bulk_update(to_update, ["gene_symbol_reported", "hgnc"])
            to_update = []

    if to_update:
        PanelAppPanelLocalCacheGeneSymbol.objects.bulk_update(to_update, ["gene_symbol_reported", "hgnc"])


class Migration(migrations.Migration):

    dependencies = [
        ("genes", "0083_one_off_canonical_gene_annotation_import_urls"),
    ]

    operations = [
        migrations.AddField(
            model_name="panelapppanellocalcachegenesymbol",
            name="hgnc",
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL,
                                    to="genes.hgnc"),
        ),
        migrations.AddField(
            model_name="panelapppanellocalcachegenesymbol",
            name="gene_symbol_reported",
            field=models.TextField(default=""),
            preserve_default=False,
        ),
        migrations.RunPython(backfill_hgnc_and_reported_symbol, migrations.RunPython.noop),
        migrations.RemoveField(
            model_name="panelapppanellocalcachegenesymbol",
            name="gene_symbol",
        ),
    ]
