from django.db import migrations, models


def set_withdrawn_record_to_other(apps, _schema_editor):
    Classification = apps.get_model('classification', 'Classification')
    for classification in Classification.objects.filter(withdrawn=True, withdraw_reason__isnull=True):
        classification.withdraw_reason = 'OTHER'
        classification.save()


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0113_uploadedclassificationsunmapped_file_type_override'),
    ]

    operations = [
        migrations.AddField(
            model_name='classification',
            name='withdraw_reason',
            field=models.CharField(blank=True, choices=[('SHARED_BY_MISTAKE', 'This record was shared by mistake'),
                                                        ('DUPLICATE', 'This record has an exact duplicate'),
                                                        ('SENSITIVE_DATA', 'This record contains sensitive data'),
                                                        ('OTHER', 'Other')], max_length=50, null=True),
        ),
        migrations.RunPython(set_withdrawn_record_to_other),
    ]
