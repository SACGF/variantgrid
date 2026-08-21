import django.db.models.deletion
from django.db import migrations, models


def _fusion_processor_genome_build(apps, _schema_editor):
    """ DRAGEN TSO500's fusion caller. Its AllFusions.csv rows are gene-level, so the VCF written from
        them has no contigs to detect a build from and the file declares none - without this the import
        stops at REQUIRES_USER_INPUT on any server with more than one build. TSO500 is run against
        hg19; when that changes, this row changes. """

    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    if grch37 := GenomeBuild.objects.filter(pk="GRCh37").first():
        VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
        VCFSourceSettings.objects.update_or_create(source_regex="^FusionProcessor",
                                                   defaults={"genome_build": grch37})


def _remove_fusion_processor_genome_build(apps, _schema_editor):
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.filter(source_regex="^FusionProcessor").delete()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0213_tag_allele_origin_bucket'),
    ]

    operations = [
        migrations.AddField(
            model_name='vcfsourcesettings',
            name='genome_build',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE,
                                    to='snpdb.genomebuild'),
        ),
        migrations.RunPython(_fusion_processor_genome_build, _remove_fusion_processor_genome_build),
    ]
