from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_vcfs_with_samples(apps):
    return apps.get_model("snpdb", "VCF").objects.filter(genotype_samples__gt=0).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0001_initial'),
        ('snpdb', '0016_api_token'),
        ('upload', '0004_file_type_manual_variant_entry_sample'),
    ]

    operations = [
        ManualOperation.operation_manage(["one_off_clear_unused_vcf_format_fields"], test=_test_has_vcfs_with_samples),
    ]
