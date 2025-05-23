# Generated by Django 4.2.9 on 2024-12-17 06:34

from django.db import migrations


def _vep_gnomad_t2t(apps, _schema_editor):
    ColumnVEPField = apps.get_model("annotation", "ColumnVEPField")
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    grch38 = GenomeBuild.objects.get(name="GRCh38")
    t2t = GenomeBuild.objects.get(name="T2T-CHM13v2.0")

    VEP_CUSTOM_GNOMAD_4 = 'o'
    VEP_CUSTOM_GNOMAD_SV = 'S'
    VEP_CUSTOM_GNOMAD_SV_NAME = 'N'

    skip_t2t = {"faf95", "faf99"}

    vep_custom = [VEP_CUSTOM_GNOMAD_4, VEP_CUSTOM_GNOMAD_SV, VEP_CUSTOM_GNOMAD_SV_NAME]
    new_records = []
    for cvf in ColumnVEPField.objects.filter(genome_build=grch38, vep_custom__in=vep_custom):
        if cvf.source_field in skip_t2t:
            continue
        cvf.pk = None
        column = cvf.column
        if column.endswith("_38"):
            column = column[:-3]
        column += "_t2t"
        cvf.column = column
        cvf.genome_build = t2t
        new_records.append(cvf)

    if new_records:
        ColumnVEPField.objects.bulk_create(new_records)


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0121_vep_columns_37_38_only"),
        ('snpdb', '0155_create_t2t_genome_build_and_contigs'),  # Where T2T GenomeBuild inserted
    ]

    operations = [
        migrations.RunPython(_vep_gnomad_t2t)
    ]
