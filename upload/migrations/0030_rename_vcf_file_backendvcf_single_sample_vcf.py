from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("upload", "0029_rename_combo_vcf_backendvcf_joint_called_vcf"),
    ]

    operations = [
        migrations.RenameField(
            model_name="backendvcf",
            old_name="vcf_file",
            new_name="single_sample_vcf",
        ),
    ]
