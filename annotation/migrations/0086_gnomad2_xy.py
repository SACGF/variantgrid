# Generated by Django 4.1.4 on 2023-12-07 11:48

from django.db import migrations

from library.django_utils import bulk_insert_class_data


def _gnomad2xy(apps, _schema_editor):
    VEP_CUSTOM_GNOMAD_2 = 'g'
    FREQUENCY_DATA = 'F'

    COLUMN_VEP_FIELD = [
        {'column': 'gnomad2_xy_ac', 'source_field_has_custom_prefix': True, 'min_vep_columns_version': 3,
         'vep_custom': VEP_CUSTOM_GNOMAD_2, 'variant_grid_column_id': 'gnomad_xy_ac', 'source_field': 'AC_male',
         'category': FREQUENCY_DATA, 'genome_build_id': 'GRCh37'},
        {'column': 'gnomad2_xy_af', 'source_field_has_custom_prefix': True, 'min_vep_columns_version': 3,
         'vep_custom': VEP_CUSTOM_GNOMAD_2, 'variant_grid_column_id': 'gnomad_xy_af', 'source_field': 'AF_male',
         'category': FREQUENCY_DATA, 'genome_build_id': 'GRCh37'},
        {'column': 'gnomad2_xy_an', 'source_field_has_custom_prefix': True, 'min_vep_columns_version': 3,
         'vep_custom': VEP_CUSTOM_GNOMAD_2, 'variant_grid_column_id': 'gnomad_xy_an', 'source_field': 'AN_male',
         'category': FREQUENCY_DATA, 'genome_build_id': 'GRCh37'},
    ]
    bulk_insert_class_data(apps, "annotation", [("ColumnVEPField", COLUMN_VEP_FIELD)])


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0085_alter_dbnsfpgeneannotationversion_options_and_more"),
    ]

    operations = [
        migrations.RunPython(_gnomad2xy)
    ]