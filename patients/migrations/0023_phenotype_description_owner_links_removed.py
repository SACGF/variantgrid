""" The PatientTextPhenotype / CohortTextPhenotype link rows go now 0022 copied them onto PhenotypeDescription,
    which holds at most one owner (#2135) """
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ("patients", "0022_phenotype_description_owners_copied"),
    ]

    operations = [
        migrations.DeleteModel(name="CohortTextPhenotype"),
        migrations.DeleteModel(name="PatientTextPhenotype"),
        migrations.AddConstraint(
            model_name="phenotypedescription",
            constraint=models.CheckConstraint(condition=models.Q(("patient__isnull", True), ("cohort__isnull", True),
                                                                 _connector="OR"),
                                              name="phenotype_description_one_owner"),
        ),
    ]
