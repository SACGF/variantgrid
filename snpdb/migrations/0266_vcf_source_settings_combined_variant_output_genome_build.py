from django.db import migrations


def _combined_variant_output_genome_build(apps, _schema_editor):
    """ The CombinedVariantOutput.tsv names no build and the splice breakpoints are positions in one, so
        without this the import fails on any server with more than one build. TSO500 is run against
        hg19, as ^FusionProcessor's row already says (0214); when that changes, both rows change. """

    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    if grch37 := GenomeBuild.objects.filter(pk="GRCh37").first():
        VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
        VCFSourceSettings.objects.update_or_create(source_regex="^DRAGEN TSO500 CombinedVariantOutput",
                                                   defaults={"genome_build": grch37})


def _remove_combined_variant_output_genome_build(apps, _schema_editor):
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.filter(source_regex="^DRAGEN TSO500 CombinedVariantOutput") \
        .update(genome_build=None)


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0265_alter_vcflengthstats_variant_class'),
    ]

    operations = [
        migrations.RunPython(_combined_variant_output_genome_build,
                             _remove_combined_variant_output_genome_build),
    ]
