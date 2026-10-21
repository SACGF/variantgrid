from django.db import migrations


def _phastcons_mammal_labels(apps, _schema_editor):
    """ Spelt like the other PhastCons / PhyloP labels ("PhyloP 46way mammal") """
    VariantGridColumn = apps.get_model("snpdb", "VariantGridColumn")
    for ways in [30, 46]:
        VariantGridColumn.objects.filter(pk=f"phastcons_{ways}_way_mammalian").update(label=f"PhastCons {ways}way mammal")


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0291_move_classification_settings'),
    ]

    operations = [
        migrations.RunPython(_phastcons_mammal_labels, migrations.RunPython.noop),
    ]
