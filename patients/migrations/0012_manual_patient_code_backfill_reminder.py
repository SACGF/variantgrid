from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_existing_patients(apps):
    Patient = apps.get_model("patients", "Patient")
    return Patient.objects.exists()


class Migration(migrations.Migration):
    dependencies = [
        ("patients", "0011_patient_patient_code"),
        ("manual", "0002_deployment"),
    ]

    operations = [
        ManualOperation.operation_other(
            args=["Patient.patient_code added (de-identified, safe to show where names were historically hidden) - if you've been using other fields as codes, manually migrate them to patient_code. See https://github.com/SACGF/variantgrid/issues/230"],
            test=_has_existing_patients,
        ),
    ]
