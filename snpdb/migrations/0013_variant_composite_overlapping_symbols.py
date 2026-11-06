from django.db import migrations


def _add_overlapping_symbols_member(apps, _schema_editor):
    """ The Variant cell names the genes an SV spans, so it needs them riding along hidden """
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    CompositeColumnMember.objects.get_or_create(composite_id="variant", column_id="overlapping_symbols",
                                                defaults={"sort_order": 9, "in_sort_menu": False})


def _remove_overlapping_symbols_member(apps, _schema_editor):
    CompositeColumnMember = apps.get_model("snpdb", "CompositeColumnMember")
    CompositeColumnMember.objects.filter(composite_id="variant", column_id="overlapping_symbols").delete()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0012_cytoband_levels'),
    ]

    operations = [
        migrations.RunPython(_add_overlapping_symbols_member, reverse_code=_remove_overlapping_symbols_member),
    ]
