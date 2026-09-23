from django.db import migrations

from manual.operations.manual_operations import ManualOperation

SOURCE_REGEX = "^DRAGEN TSO500 CombinedVariantOutput"


def _remove_combined_variant_output_source_settings(apps, _schema_editor):
    """ The CombinedVariantOutput import writes no VCF any more (#1904), so nothing is called against
        this source: it bound the read-support FORMAT fields of the record-less VCF (0264) and the
        build it declared nowhere (0266) """
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.filter(source_regex=SOURCE_REGEX).delete()


def _restore_combined_variant_output_source_settings(apps, _schema_editor):
    GenomeBuild = apps.get_model("snpdb", "GenomeBuild")
    VCFSourceSettings = apps.get_model("snpdb", "VCFSourceSettings")
    VCFSourceSettings.objects.update_or_create(
        source_regex=SOURCE_REGEX,
        defaults={
            "genome_build": GenomeBuild.objects.filter(pk="GRCh37").first(),
            "sample_field_overrides": {
                "allele_depth_field": None,
                "alt_depth_field": "ALT_READS",
                "ref_depth_field": "REF_READS",
                "read_depth_field": None,
                "allele_frequency_field": None,
            },
        })


def _has_record_less_combined_variant_output_vcfs(apps):
    VCF = apps.get_model("snpdb", "VCF")
    return VCF.objects.filter(source__startswith="DRAGEN TSO500 CombinedVariantOutput").exists()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0274_tag_config_show_resolved_variant_tags'),
    ]

    operations = [
        migrations.RunPython(_remove_combined_variant_output_source_settings,
                             _restore_combined_variant_output_source_settings),
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_delete_combined_variant_output_vcfs"]),
                        note="Soft-delete the record-less VCFs earlier CombinedVariantOutput imports wrote, "
                             "so they leave the data and run pages (#1904)",
                        test=_has_record_less_combined_variant_output_vcfs),
    ]
