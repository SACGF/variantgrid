# Generated by Django 4.2.5 on 2023-10-31 04:57

from django.db import migrations, models


def _set_allow_custom_strengths(apps, schema):
    EvidenceKey = apps.get_model('classification', 'EvidenceKey')
    EvidenceKey.objects.filter(value_type="C").update(crit_allows_override_strengths=True)

def _dummy_reverse(apps, schema):
    pass

class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0116_merge_20231027_1207'),
    ]

    operations = [
        migrations.AddField(
            model_name='evidencekey',
            name='crit_allows_override_strengths',
            field=models.BooleanField(blank=True, default=False),
        ),
        migrations.AddField(
            model_name='evidencekey',
            name='crit_uses_points',
            field=models.BooleanField(blank=True, default=False),
        ),
        migrations.RunPython(_set_allow_custom_strengths, _dummy_reverse)
    ]