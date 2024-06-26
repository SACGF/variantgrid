# Generated by Django 4.2.2 on 2023-12-20 00:14
from django.db import migrations

def _one_off_annotation_version_use_cosmic_vcf(apps, _schema_editor):
    """ We have long used Cosmic VCFs, but used the VEP cosmic version in VariantAnnotationVersion
        so we need to update historical versions
     """

    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    VariantAnnotationVersion.objects.filter(columns_version=1).update(cosmic=95)
    VariantAnnotationVersion.objects.filter(columns_version=2).update(cosmic=97)
    VariantAnnotationVersion.objects.filter(columns_version=3).update(cosmic=99)


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0086_gnomad2_xy'),
    ]

    operations = [
        migrations.RunPython(_one_off_annotation_version_use_cosmic_vcf)
    ]
