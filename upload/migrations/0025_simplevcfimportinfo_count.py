# Generated by Django 4.2.18 on 2025-02-04 03:48

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('upload', '0024_modifiedimportedvariant_operation'),
    ]

    operations = [
        migrations.AddField(
            model_name='simplevcfimportinfo',
            name='count',
            field=models.IntegerField(default=0),
        ),
    ]
