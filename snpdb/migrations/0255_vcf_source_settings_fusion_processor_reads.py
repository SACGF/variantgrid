from django.db import migrations


def _fusion_processor_read_support(apps, _schema_editor):
    """ The VCF written from AllFusions.csv carries the sample's read support rather than a genotype:
        ALT_READS are the reads supporting the fusion and REF_READS those across the junctions that
        don't. Bound as alt and ref depth so the sample node's min-reads threshold and VAF (derived as
        alt / (alt + ref)) work, the way ^SpliceGirl's AD and DP are. VCFs imported while the writer
        still declared a GT of ./. have no genotype either, so the field is cleared on them. """
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.update_or_create(
        source_regex="^FusionProcessor",
        defaults={
            "sample_field_overrides": {
                "allele_depth_field": None,
                "alt_depth_field": "ALT_READS",  # reads supporting the fusion
                "ref_depth_field": "REF_READS",  # reads across the junctions that do not
                "read_depth_field": None,
                "allele_frequency_field": None,
            },
        },
    )
    VCF = apps.get_model("snpdb", "VCF")
    VCF.objects.filter(source__startswith="FusionProcessor").update(genotype_field=None)


def _remove_fusion_processor_read_support(apps, _schema_editor):
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.filter(source_regex="^FusionProcessor").update(sample_field_overrides={})


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0254_one_off_backfill_vcf_genotype_field'),
    ]

    operations = [
        migrations.RunPython(_fusion_processor_read_support, _remove_fusion_processor_read_support),
    ]
