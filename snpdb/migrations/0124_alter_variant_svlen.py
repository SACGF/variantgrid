# Generated by Django 4.2.10 on 2024-04-04 00:46

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0123_new_vg_columns_and_custom_columns'),
    ]

    operations = [
        migrations.AlterField(
            model_name='variant',
            name='svlen',
            field=models.IntegerField(blank=True, null=True),
        ),
    ]
