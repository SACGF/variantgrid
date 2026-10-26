from django.db import migrations


class Migration(migrations.Migration):
    """ BackendVCF and the DRAGEN TSO500 UploadData move to seqauto (#2490), which keeps their tables -
        state only, seqauto/migrations/0005_models_from_upload.py takes them over """

    dependencies = [
        ('seqauto', '0004_seed_data'),
        ('upload', '0002_modified_imported_variant_conflicting_duplicate'),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.DeleteModel(name='BackendVCF'),
                migrations.DeleteModel(name='UploadedDragenTSO500CombinedVariantOutput'),
                migrations.DeleteModel(name='UploadedDragenTSO500MetricsOutput'),
            ],
        ),
    ]
