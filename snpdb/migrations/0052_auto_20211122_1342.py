# Generated by Django 3.2.6 on 2021-11-22 03:12

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0051_auto_20211122_1342'),
    ]

    operations = [
        migrations.AddField(
            model_name='lab',
            name='country',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.PROTECT, to='snpdb.country'),
        ),
        migrations.AddField(
            model_name='lab',
            name='state',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.PROTECT, to='snpdb.state'),
        ),
    ]
