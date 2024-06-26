# Generated by Django 3.1.3 on 2021-03-09 04:32

from django.db import migrations

from library.django_utils import bulk_insert_class_data


def _new_vep_annotation_spliceai_gene(apps, _schema_editor):
    COLUMN_VEP_FIELD = [
        {'column': 'spliceai_gene_symbol', 'vep_plugin': 'a', 'source_field_has_custom_prefix': False,
         'source_field_processing_description': None, 'vep_custom': None,
         'variant_grid_column_id': 'spliceai_gene_symbol', 'source_field': 'SpliceAI_pred_SYMBOL', 'category': 'S'},
    ]
    bulk_insert_class_data(apps, "annotation", [("ColumnVEPField", COLUMN_VEP_FIELD)])


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0026_auto_20210303_1554'),
        ('snpdb', '0023_new_column_spliaceai_gene'),
    ]

    operations = [
        migrations.RunPython(_new_vep_annotation_spliceai_gene)
    ]
