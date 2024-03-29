# Generated by Django 3.2.1 on 2021-05-19 07:40

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0033_one_off_fix_drug_response_typo'),
    ]

    operations = [
        migrations.AlterField(
            model_name='nodecountsettings',
            name='built_in_filter',
            field=models.CharField(choices=[('T', 'Total'), ('C', 'ClinVar LP/P'), ('O', 'OMIM Phenotype'), ('I', 'High or Mod impact'), ('G', 'Classified'), ('P', 'Classified Pathogenic'), ('M', 'COSMIC')], max_length=1, null=True),
        ),
    ]
