from django.db import migrations


def _delete_refseq_gene_cached_web_resources(apps, _schema_editor):
    """ RefSeq gene summaries and NCBI gene symbol aliases come from the cdot import (#2431) """
    CachedWebResource = apps.get_model("snpdb", "CachedWebResource")
    CachedWebResource.objects.filter(name__in=["RefSeq Gene Info", "RefSeq Gene Summary"]).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("genes", "0003_seed_data"),
        ("snpdb", "0003_utr_annotator_avi_splicevault_columns"),
    ]

    operations = [
        migrations.RunPython(_delete_refseq_gene_cached_web_resources, migrations.RunPython.noop),
    ]
