# Generated by Django 4.2.10 on 2024-05-10 06:19

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0071_rename_oe_lof_gnomadgeneconstraint_lof_oe_and_more'),
        ('analysis', '0081_one_off_modify_sample_tab_auto_analysis_template'),
    ]

    operations = [
        migrations.AddField(
            model_name='karyomappinggene',
            name='gene_symbol',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='genes.genesymbol'),
        ),
    ]
