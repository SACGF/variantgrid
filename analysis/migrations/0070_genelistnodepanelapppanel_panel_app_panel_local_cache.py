# Generated by Django 4.1.3 on 2023-01-04 00:48

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0063_panelapppanellocalcache_and_more'),
        ('analysis', '0069_remove_genelistnodepanelapppanel_panel_app_panel_local_cache_gene_list_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='genelistnodepanelapppanel',
            name='panel_app_panel_local_cache',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='genes.panelapppanellocalcache'),
        ),
    ]
