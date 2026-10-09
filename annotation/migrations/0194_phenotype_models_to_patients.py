""" The phenotype text models moved to patients (#2135) - patients 0020 renamed their tables, so they only leave
    annotation's state here """
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0193_one_off_match_patient_phenotypes_stale_fuzzy"),
        ("patients", "0020_phenotype_models_from_annotation"),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.DeleteModel(name="CohortTextPhenotype"),
                migrations.DeleteModel(name="PatientTextPhenotype"),
                migrations.DeleteModel(name="TextPhenotypeMatch"),
                migrations.DeleteModel(name="TextPhenotypeSentence"),
                migrations.DeleteModel(name="TextPhenotype"),
                migrations.DeleteModel(name="PhenotypeMatchVersion"),
                migrations.DeleteModel(name="PhenotypeDescription"),
                migrations.DeleteModel(name="DescriptionProcessingStatus"),
            ],
        ),
    ]
