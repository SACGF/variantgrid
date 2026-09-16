from django.db import migrations

# SACGF/variantgrid#1835 - a splice label is now canonical: lower-case tokens joined by underscores,
# the one form every way of writing a junction's name resolves to (@see genes.gene_splice).
# 0093 seeded the labels as a report writes them, so the rows it made are re-pointed here; the
# display, which is what a report gets, is unchanged.
CANONICAL_LABELS = {
    ("AR", "V7"): "v_7",
    ("EGFR", "vIII"): "v_iii",
    ("MET", "ex14skip"): "exon_14_skipping",
}


def _canonicalise_labels(apps, _schema_editor):
    SpliceEvent = apps.get_model("genes", "SpliceEvent")
    for (gene_symbol, label), canonical in CANONICAL_LABELS.items():
        SpliceEvent.objects.filter(gene_symbol=gene_symbol, label=label).update(label=canonical)


def _restore_seeded_labels(apps, _schema_editor):
    SpliceEvent = apps.get_model("genes", "SpliceEvent")
    for (gene_symbol, label), canonical in CANONICAL_LABELS.items():
        SpliceEvent.objects.filter(gene_symbol=gene_symbol, label=canonical).update(label=label)


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0094_alter_spliceevent_label'),
    ]

    operations = [
        migrations.RunPython(_canonicalise_labels, reverse_code=_restore_seeded_labels),
    ]
