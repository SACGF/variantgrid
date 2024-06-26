# Generated by Django 4.2.9 on 2024-01-22 06:24

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0077_alter_populationnodegnomadpopulation_population'),
    ]

    operations = [
        migrations.AddField(
            model_name='damagenode',
            name='alphamissense_rankscore_allow_null',
            field=models.BooleanField(default=True),
        ),
        migrations.AddField(
            model_name='damagenode',
            name='alphamissense_rankscore_min',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='damagenode',
            name='alphamissense_rankscore_required',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='populationnode',
            name='gnomad_fafmax_faf95_max',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='populationnode',
            name='gnomad_fafmax_faf99_max',
            field=models.BooleanField(default=False),
        ),
    ]
