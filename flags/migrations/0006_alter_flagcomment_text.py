# Generated by Django 4.0.4 on 2022-04-29 01:01

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('flags', '0005_no_public_sharing_flag'),
    ]

    operations = [
        migrations.AlterField(
            model_name='flagcomment',
            name='text',
            field=models.TextField(blank=True),
        ),
    ]
