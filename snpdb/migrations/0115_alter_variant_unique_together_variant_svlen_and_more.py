# Generated by Django 4.2.9 on 2024-02-14 03:10

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0114_alter_variant_end'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='variant',
            unique_together=set(),
        ),
        migrations.AddField(
            model_name='variant',
            name='svlen',
            field=models.IntegerField(null=True),
        ),
        migrations.AlterUniqueTogether(
            name='variant',
            unique_together={('locus', 'alt', 'svlen')},
        ),
    ]