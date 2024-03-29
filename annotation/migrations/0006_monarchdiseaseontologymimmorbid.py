# Generated by Django 3.1 on 2020-11-06 04:02

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0005_monarchdiseaseontologygenerelationship'),
    ]

    operations = [
        migrations.CreateModel(
            name='MonarchDiseaseOntologyMIMMorbid',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('relationship', models.TextField()),
                ('omim_id', models.IntegerField()),
                ('mondo', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='annotation.monarchdiseaseontology')),
            ],
            options={
                'unique_together': {('mondo', 'relationship', 'omim_id')},
            },
        ),
    ]
