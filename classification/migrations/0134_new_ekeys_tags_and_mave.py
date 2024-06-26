# Generated by Django 4.2.10 on 2024-04-04 01:52

from django.db import migrations


def _new_ekey_tags_and_hgvs_g(apps, _schema_editor):
    EvidenceKey = apps.get_model("classification", "EvidenceKey")

    # This is done in annotation now
    EvidenceKey.objects.filter(key='g_hgvs').update(variantgrid_column_id='hgvs_g')

    # New one for tags
    OTHER_DATABASE = 'DB'
    FREE_ENTRY = 'F'

    kwargs = {
        'key': 'variant_tags', 'mandatory': False, 'max_share_level': 'logged_in_users', 'order': 2,
        'label': 'Variant Tags', 'sub_label': None, 'description': 'Variant Tags',
        'examples': ['Artifact x 2'], 'options': [], 'see': None, 'evidence_category': OTHER_DATABASE,
        'value_type': FREE_ENTRY, 'default_crit_evaluation': None, 'allow_custom_values': False, 'hide': True,
        'immutable': False, 'copy_consensus': False, 'variantgrid_column_id': None,
    }
    EvidenceKey.objects.create(**kwargs)


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0133_fix_om3_om4'),
        ('snpdb', '0123_new_vg_columns_and_custom_columns'),
    ]

    operations = [
        migrations.RunPython(_new_ekey_tags_and_hgvs_g)
    ]
