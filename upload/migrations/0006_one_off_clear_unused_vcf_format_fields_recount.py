from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_vcfs_with_samples(apps):
    return apps.get_model("snpdb", "VCF").objects.filter(genotype_samples__gt=0).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('upload', '0005_one_off_clear_unused_vcf_format_fields'),
    ]

    operations = [
        # Registered again: a deployment that already ran it unbound the fields without recounting sample stats
        ManualOperation.operation_manage(["one_off_clear_unused_vcf_format_fields"],
                                         note="Recounts sample stats for VCFs whose unused GT field was unbound, so "
                                              "their contents show (#2833)",
                                         test=_test_has_vcfs_with_samples),
    ]
