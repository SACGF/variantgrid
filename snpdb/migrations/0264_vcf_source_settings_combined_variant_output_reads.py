from django.db import migrations


def _combined_variant_output_read_support(apps, _schema_editor):
    """ The VCF written from a TSO 500 CombinedVariantOutput's '[Splice Variants]' carries the
        sample's read support rather than a genotype: ALT_READS are the reads crossing the junction
        and REF_READS those across the reference transcript. Bound as alt and ref depth so the sample
        node's min-reads threshold and VAF (derived as alt / (alt + ref)) work, the way ^FusionProcessor's
        are - the frequency it yields is the junction ratio. """
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.update_or_create(
        source_regex="^DRAGEN TSO500 CombinedVariantOutput",
        defaults={
            "sample_field_overrides": {
                "allele_depth_field": None,
                "alt_depth_field": "ALT_READS",  # reads supporting the splice junction
                "ref_depth_field": "REF_READS",  # reads across the reference transcript
                "read_depth_field": None,
                "allele_frequency_field": None,
            },
        },
    )


def _remove_combined_variant_output_read_support(apps, _schema_editor):
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.filter(source_regex="^DRAGEN TSO500 CombinedVariantOutput") \
        .update(sample_field_overrides={})


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0263_mandatory_columns_all_collections'),
    ]

    operations = [
        migrations.RunPython(_combined_variant_output_read_support,
                             _remove_combined_variant_output_read_support),
    ]
