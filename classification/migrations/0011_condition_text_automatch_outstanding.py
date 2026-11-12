from django.db import migrations


def _flag_outstanding_for_automatch(apps, _schema_editor):
    """ Automatch now accepts exact MONDO synonyms, OMIM titles converted to MONDO and a trailing mode of inheritance
        (SACGF/variantgrid_private#3948): queue texts with unmatched classifications for the beat sweep to retry """
    ConditionText = apps.get_model("classification", "ConditionText")
    ConditionText.objects.filter(classifications_count_outstanding__gt=0).update(pending_automatch=True)


class Migration(migrations.Migration):

    dependencies = [
        ("classification", "0010_merge_20261102_1217"),
    ]

    operations = [
        migrations.RunPython(_flag_outstanding_for_automatch, migrations.RunPython.noop),
    ]
