import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('patients', '0013_extraction_specimen_surrogate_pk'),
        ('seqauto', '0045_jointcalledvcf_sequencing_samples'),
    ]

    operations = [
        migrations.AddField(
            model_name='sequencingsample',
            name='extraction',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL,
                                    to='patients.extraction'),
        ),
    ]
