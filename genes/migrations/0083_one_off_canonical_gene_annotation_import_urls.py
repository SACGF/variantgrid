from django.db import migrations

from genes.gene_annotation_import_urls import canonical_import_url


def _canonicalise_import_urls(apps, _schema_editor):
    """ GeneAnnotationImport.url is how we recognise an already imported GTF/GFF, so the same file
        recorded under a different scheme (ftp/http vs https) or format (Ensembl gff3 vs gtf) created a
        second row. Collapse each set onto one row so URL comparison identifies a release. """
    GeneAnnotationImport = apps.get_model("genes", "GeneAnnotationImport")
    GeneVersion = apps.get_model("genes", "GeneVersion")
    TranscriptVersion = apps.get_model("genes", "TranscriptVersion")
    GeneAnnotationRelease = apps.get_model("genes", "GeneAnnotationRelease")

    by_canonical_url = {}
    for gai in GeneAnnotationImport.objects.all():
        by_canonical_url.setdefault(canonical_import_url(gai.url), []).append(gai)

    for canonical_url, imports in by_canonical_url.items():
        # Prefer a row that already has the canonical URL, so we don't renumber needlessly
        already_canonical = [gai for gai in imports if gai.url == canonical_url]
        keeper = min(already_canonical or imports, key=lambda gai: gai.pk)

        for gai in imports:
            if gai.pk == keeper.pk:
                continue
            GeneVersion.objects.filter(import_source=gai).update(import_source=keeper)
            TranscriptVersion.objects.filter(import_source=gai).update(import_source=keeper)
            GeneAnnotationRelease.objects.filter(gene_annotation_import=gai).update(gene_annotation_import=keeper)
            gai.delete()

        if keeper.url != canonical_url:
            keeper.url = canonical_url
            keeper.save()


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0082_one_off_panel_app_deleted_status_to_flag'),
    ]

    operations = [
        migrations.RunPython(_canonicalise_import_urls, migrations.RunPython.noop),
    ]
