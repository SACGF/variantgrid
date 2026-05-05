from django.db import migrations


def _reminder(apps, schema_editor):
    """
    MANUAL MIGRATION REQUIRED for sites that were using Patient.last_name as a
    de-identified code (issue #230 workaround).

    The new Patient.patient_code field is the canonical de-identified display
    label. It is intentionally NOT auto-populated from last_name, because
    last_name may legitimately contain real surnames on some deployments and a
    blind copy would leak PII into a field designed to be shown everywhere.

    Each deployment must decide and run its own backfill, e.g.:

        # If you know last_name on this deployment is always a code:
        UPDATE patients_patient
           SET patient_code = last_name
         WHERE patient_code IS NULL
           AND last_name IS NOT NULL;

        # If only some labs/imports used last_name as a code, scope the UPDATE
        # by the relevant Sample/VCF/Lab join.

    Until the backfill runs, %(patient_code)s in grid_sample_label_template
    will resolve to '' for legacy patients.
    """


class Migration(migrations.Migration):
    dependencies = [("patients", "0011_patient_patient_code")]
    operations = [migrations.RunPython(_reminder, migrations.RunPython.noop)]
