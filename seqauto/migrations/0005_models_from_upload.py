import django.db.models.deletion
from django.db import migrations, models

MOVED_MODELS = ['backendvcf', 'uploadeddragentso500combinedvariantoutput', 'uploadeddragentso500metricsoutput']


def _move_content_types(apps, schema_editor):
    """ Keeps the moved models' permissions - a fresh database has no content types yet """
    ContentType = apps.get_model('contenttypes', 'ContentType')
    existing = set(ContentType.objects.filter(app_label='seqauto', model__in=MOVED_MODELS).values_list('model', flat=True))
    ContentType.objects.filter(app_label='upload', model__in=set(MOVED_MODELS) - existing).update(app_label='seqauto')


class Migration(migrations.Migration):
    """ Takes over the tables of upload's BackendVCF and DRAGEN TSO500 UploadData (#2490) - state only """

    dependencies = [
        ('contenttypes', '0002_remove_content_type_name'),
        ('seqauto', '0004_seed_data'),
        ('upload', '0003_move_sequencing_models_to_seqauto'),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.CreateModel(
                    name='BackendVCF',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('joint_called_vcf', models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE, to='seqauto.jointcalledvcf')),
                        ('single_sample_vcf', models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE, to='seqauto.singlesamplevcf')),
                        ('uploaded_vcf', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='upload.uploadedvcf')),
                    ],
                    options={
                        'db_table': 'upload_backendvcf',
                    },
                ),
                migrations.CreateModel(
                    name='UploadedDragenTSO500CombinedVariantOutput',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('file_upload', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='upload.fileupload')),
                    ],
                    options={
                        'db_table': 'upload_uploadeddragentso500combinedvariantoutput',
                    },
                ),
                migrations.CreateModel(
                    name='UploadedDragenTSO500MetricsOutput',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('file_upload', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='upload.fileupload')),
                    ],
                    options={
                        'db_table': 'upload_uploadeddragentso500metricsoutput',
                    },
                ),
            ],
        ),
        migrations.RunPython(_move_content_types, migrations.RunPython.noop),
    ]
