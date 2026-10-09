""" A PhenotypeDescription points at the patient or cohort that owns it, so deleting the owner deletes it (#2135).
    The owner and approval come off the PatientTextPhenotype / CohortTextPhenotype link rows, which go; a
    description nothing holds (left by deleted patients) is deleted with its sentences. """
import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models

COPY_OWNERS_SQL = [
    """UPDATE patients_phenotypedescription d
       SET patient_id = l.patient_id, approved_by_id = l.approved_by_id
       FROM patients_patienttextphenotype l
       WHERE l.phenotype_description_id = d.id""",
    """UPDATE patients_phenotypedescription d
       SET cohort_id = l.cohort_id, approved_by_id = l.approved_by_id
       FROM patients_cohorttextphenotype l
       WHERE l.phenotype_description_id = d.id""",
    # Check the deferred FKs now, or the ALTER TABLEs below fail on pending trigger events
    "SET CONSTRAINTS ALL IMMEDIATE",
]

COPY_OWNERS_REVERSE_SQL = [
    """INSERT INTO patients_patienttextphenotype (patient_id, phenotype_description_id, approved_by_id)
       SELECT patient_id, id, approved_by_id FROM patients_phenotypedescription WHERE patient_id IS NOT NULL""",
    """INSERT INTO patients_cohorttextphenotype (cohort_id, phenotype_description_id, approved_by_id)
       SELECT cohort_id, id, approved_by_id FROM patients_phenotypedescription WHERE cohort_id IS NOT NULL""",
]


def _delete_unowned_descriptions(apps, schema_editor):
    """ Unowned = no patient, no cohort and nothing else pointing at it - any other app's relation in the
        migration state (an SA Path request's link rows) holds a description too """
    PhenotypeDescription = apps.get_model("patients", "PhenotypeDescription")
    TextPhenotypeSentence = apps.get_model("patients", "TextPhenotypeSentence")
    quote_name = schema_editor.connection.ops.quote_name

    held_conditions = ["d.patient_id IS NOT NULL", "d.cohort_id IS NOT NULL"]
    for relation in PhenotypeDescription._meta.related_objects:
        if relation.related_model is TextPhenotypeSentence:
            continue
        table = quote_name(relation.related_model._meta.db_table)
        column = quote_name(relation.field.column)
        held_conditions.append(f"EXISTS (SELECT 1 FROM {table} r WHERE r.{column} = d.id)")
    unowned_sql = f"SELECT d.id FROM patients_phenotypedescription d WHERE NOT ({' OR '.join(held_conditions)})"

    with schema_editor.connection.cursor() as cursor:
        cursor.execute(f"DELETE FROM patients_textphenotypesentence WHERE phenotype_description_id IN ({unowned_sql})")
        cursor.execute(f"DELETE FROM patients_phenotypedescription WHERE id IN ({unowned_sql})")
        cursor.execute("SET CONSTRAINTS ALL IMMEDIATE")


class Migration(migrations.Migration):

    dependencies = [
        # After annotation's state lets go of the models, so any other app's links to them are in the state
        ("annotation", "0193_phenotype_models_to_patients"),
        ("patients", "0020_phenotype_models_from_annotation"),
        ("snpdb", "0281_genomic_intervals_collection_error_message"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.AddField(
            model_name="phenotypedescription",
            name="patient",
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE,
                                       related_name="phenotype_description", to="patients.patient"),
        ),
        migrations.AddField(
            model_name="phenotypedescription",
            name="cohort",
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE,
                                       related_name="phenotype_description", to="snpdb.cohort"),
        ),
        migrations.AddField(
            model_name="phenotypedescription",
            name="approved_by",
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL,
                                    to=settings.AUTH_USER_MODEL),
        ),
        migrations.RunSQL(COPY_OWNERS_SQL, reverse_sql=COPY_OWNERS_REVERSE_SQL),
        migrations.DeleteModel(name="CohortTextPhenotype"),
        migrations.DeleteModel(name="PatientTextPhenotype"),
        migrations.RunPython(_delete_unowned_descriptions, migrations.RunPython.noop),
        migrations.AddConstraint(
            model_name="phenotypedescription",
            constraint=models.CheckConstraint(condition=models.Q(("patient__isnull", True), ("cohort__isnull", True),
                                                                 _connector="OR"),
                                              name="phenotype_description_one_owner"),
        ),
    ]
