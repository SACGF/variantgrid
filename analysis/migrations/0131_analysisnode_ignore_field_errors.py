from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0130_samplenodesamplefilter'),
    ]

    operations = [
        migrations.AddField(
            model_name='analysisnode',
            name='ignore_field_errors',
            field=models.BooleanField(default=False),
        ),
    ]
