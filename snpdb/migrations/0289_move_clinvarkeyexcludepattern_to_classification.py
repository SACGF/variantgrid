# State-only: the model moves to classification, which keeps using the table snpdb_clinvarkeyexcludepattern (#2276)

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0288_one_off_allele_liftover_stuck_processing'),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.DeleteModel(
                    name='ClinVarKeyExcludePattern',
                ),
            ],
            database_operations=[],
        ),
    ]
