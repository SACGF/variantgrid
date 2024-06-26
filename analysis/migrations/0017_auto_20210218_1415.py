# Generated by Django 3.1.3 on 2021-02-18 03:45

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('pathtests', '0004_auto_20210121_1754'),
        ('analysis', '0016_one_off_fix_liftover_existing_variant_tags'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='genelistnode',
            name='pathology_test_gene_list',
        ),
        migrations.AddField(
            model_name='genelistnode',
            name='pathology_test_version',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='pathtests.pathologytestversion'),
        ),
    ]
