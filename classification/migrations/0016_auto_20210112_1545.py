# Generated by Django 3.1 on 2021-01-12 05:15

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0015_auto_20210108_1435'),
    ]

    operations = [
        migrations.DeleteModel(
            name='ConditionAliasSearchCache',
        ),
        migrations.DeleteModel(
            name='ConditionAlias',
        )
    ]