from django.db import migrations, models

# SACGF/variantgrid#1909 - the CIViC variant for each named junction CIViC records. EGFRvII has none.
CIVIC_VARIANT_IDS = {
    ("MET", "exon_14_skipping"): 324,
    ("AR", "v_7"): 362,
    ("EGFR", "v_iii"): 312,
}


def _seed_civic_variant_ids(apps, _schema_editor):
    SpliceEvent = apps.get_model("genes", "SpliceEvent")
    for (gene_symbol, label), civic_variant_id in CIVIC_VARIANT_IDS.items():
        SpliceEvent.objects.filter(gene_symbol=gene_symbol, label=label).update(civic_variant_id=civic_variant_id)


def _clear_civic_variant_ids(apps, _schema_editor):
    SpliceEvent = apps.get_model("genes", "SpliceEvent")
    SpliceEvent.objects.update(civic_variant_id=None)


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0101_seed_egfr_vii_splice_event'),
    ]

    operations = [
        migrations.AddField(
            model_name='spliceevent',
            name='civic_variant_id',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.RunPython(_seed_civic_variant_ids, reverse_code=_clear_civic_variant_ids),
    ]
