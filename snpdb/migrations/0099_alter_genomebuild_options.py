# Generated by Django 4.2.2 on 2023-09-05 03:54

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0098_columns_collections_add_internally_classified_labs'),
    ]

    operations = [
        migrations.AlterModelOptions(
            name='genomebuild',
            options={'base_manager_name': 'objects', 'ordering': ['name']},
        ),
    ]
