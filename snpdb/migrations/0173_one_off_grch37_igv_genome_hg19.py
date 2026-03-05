from django.db import migrations


def set_grch37_igv_genome_hg19(apps, schema_editor):
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    GenomeBuild.objects.filter(name="GRCh37", igv_genome="b37").update(igv_genome="hg19")


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0172_alter_sample_variants_type_and_more'),
    ]

    operations = [
        migrations.RunPython(set_grch37_igv_genome_hg19, migrations.RunPython.noop),
    ]
