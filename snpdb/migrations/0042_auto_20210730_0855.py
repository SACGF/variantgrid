# Generated by Django 3.1.6 on 2021-07-29 23:25

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0041_merge_20210730_0851'),
    ]

    operations = [
        migrations.AddField(
            model_name='clinvarkey',
            name='assertion_method_lookup',
            field=models.JSONField(default=dict),
        ),
        migrations.AddField(
            model_name='clinvarkey',
            name='default_affected_status',
            field=models.TextField(blank=True, choices=[('yes', 'Yes'), ('no', 'No'), ('unknown', 'Unknown'), ('not provided', 'Not Provided'), ('not applicable', 'Not Applicable')], null=True),
        ),
        migrations.DeleteModel(
            name='ClinVarKeyAssertionMethod',
        ),
    ]