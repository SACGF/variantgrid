# Generated by Django 3.1.3 on 2021-03-24 04:25

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('seqauto', '0010_auto_20210324_1455'),
        ('upload', '0005_simplevcfimportinfo'),
    ]

    operations = [
        migrations.AlterField(
            model_name='backendvcf',
            name='combo_vcf',
            field=models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE, to='seqauto.samplesheetcombinedvcffile2'),
        ),
        migrations.AlterField(
            model_name='backendvcf',
            name='vcf_file',
            field=models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE, to='seqauto.vcffile2'),
        ),
    ]