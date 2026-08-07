from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0037_backfill_vcffromsequencingrun_variant_caller"),
        ("upload", "0030_rename_vcf_file_backendvcf_single_sample_vcf"),
    ]

    operations = [
        migrations.RenameModel(
            old_name="VCFFile",
            new_name="SingleSampleVCF",
        ),
    ]
