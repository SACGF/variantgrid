# Generated by Django 3.1.3 on 2021-06-22 05:29

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0035_merge_20210520_1246'),
        ('annotation', '0033_one_off_impact_change_order'),
    ]

    operations = [
        migrations.AlterField(
            model_name='clinvar',
            name='variant',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='snpdb.variant'),
        ),
    ]
