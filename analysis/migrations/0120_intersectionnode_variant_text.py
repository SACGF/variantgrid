import django.contrib.postgres.fields
from django.db import migrations, models

# Panel indices before this migration - the HGVS and custom range panels merge into one free-text
# panel, which renumbers the panels after them
OLD_SELECTED_INTERVALS = 0
OLD_CUSTOM_INTERVAL = 1
OLD_HGVS = 2
OLD_BACKEND_ENRICHMENT_KIT = 3
OLD_CONTIG = 4

NEW_PANEL = {
    OLD_SELECTED_INTERVALS: 0,
    OLD_CUSTOM_INTERVAL: 1,
    OLD_HGVS: 1,
    OLD_BACKEND_ENRICHMENT_KIT: 2,
    OLD_CONTIG: 3,
}


def _convert_to_variant_text(apps, _schema_editor):
    IntersectionNode = apps.get_model("analysis", "IntersectionNode")
    GenomicInterval = apps.get_model("snpdb", "GenomicInterval")

    genomic_interval_ids = []
    records = []
    qs = IntersectionNode.objects.all().select_related("genomic_interval")
    for node in qs.iterator():
        lines = []
        if genomic_interval := node.genomic_interval:
            genomic_interval_ids.append(genomic_interval.pk)
            region = f"{genomic_interval.chrom}:{genomic_interval.start}-{genomic_interval.end}"
            lines.append(region)
            node.variant_regions = [region]
        if hgvs_string := node.hgvs_string:
            lines.append(hgvs_string)
            if node.hgvs_variant_id:
                node.variant_ids = [node.hgvs_variant_id]
            else:
                node.variant_text_unresolved = [hgvs_string]
        if lines:
            node.variant_text = "\n".join(lines)
        node.accordion_panel = NEW_PANEL.get(node.accordion_panel, node.accordion_panel)
        records.append(node)

    if records:
        IntersectionNode.objects.bulk_update(
            records, ["variant_text", "variant_ids", "variant_regions", "variant_text_unresolved", "accordion_panel"],
            batch_size=1000)
    # These were created by the custom range panel and nothing else points at them
    GenomicInterval.objects.filter(pk__in=genomic_interval_ids).delete()


class Migration(migrations.Migration):

    dependencies = [
        ("analysis", "0119_damagenode_variant_class"),
    ]

    operations = [
        migrations.AddField(
            model_name="intersectionnode",
            name="variant_ids",
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.IntegerField(), blank=True, default=list
            ),
        ),
        migrations.AddField(
            model_name="intersectionnode",
            name="variant_regions",
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.TextField(), blank=True, default=list
            ),
        ),
        migrations.AddField(
            model_name="intersectionnode",
            name="variant_text",
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name="intersectionnode",
            name="variant_text_unresolved",
            field=django.contrib.postgres.fields.ArrayField(
                base_field=models.TextField(), blank=True, default=list
            ),
        ),
        migrations.RunPython(_convert_to_variant_text, migrations.RunPython.noop),
    ]
