"""
Extraction gets a local reference beside external_pk, which is nullable because not every deployment
has an external system managing it - the same reason Patient carries patient_code and Specimen its own
reference_id.
"""

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('patients', '0014_extraction_date_and_timestamps'),
    ]

    operations = [
        migrations.AddField(
            model_name='extraction',
            name='reference_id',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AlterUniqueTogether(
            name='extraction',
            unique_together={('specimen', 'reference_id')},
        ),
    ]
