# Generated by Django 4.2.11 on 2024-06-24 01:01

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0141_populate_somatic_sort_order'),
    ]

    operations = [
        migrations.AddField(
            model_name='importedalleleinfo',
            name='dirty_message',
            field=models.TextField(blank=True, null=True),
        ),
    ]
