from django.db import migrations


def _sample_patient_from_extraction(apps, schema_editor):
    """ Samples linked to an extraction before Sample.save filled in its patient """
    with schema_editor.connection.cursor() as cursor:
        cursor.execute("""
            UPDATE snpdb_sample s
            SET patient_id = sp.patient_id
            FROM patients_extraction e
            JOIN patients_specimen sp ON sp.id = e.specimen_id
            WHERE s.extraction_id = e.id
              AND s.patient_id IS NULL
        """)


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0266_vcf_source_settings_combined_variant_output_genome_build"),
        ("patients", "0018_alter_patient_last_name"),
    ]

    operations = [
        migrations.RunPython(_sample_patient_from_extraction, reverse_code=migrations.RunPython.noop),
    ]
