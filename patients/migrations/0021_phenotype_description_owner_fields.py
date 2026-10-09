""" A PhenotypeDescription points at the patient or cohort that owns it, so deleting the owner deletes it (#2135).
    0022 copies the owners and approval off the PatientTextPhenotype / CohortTextPhenotype link rows, 0023 drops them. """
import django.db.models.deletion
from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        # After annotation's state lets go of the models, so any other app's links to them are in the state
        ("annotation", "0194_phenotype_models_to_patients"),
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
    ]
