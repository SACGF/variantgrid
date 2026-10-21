from django.db import migrations


class Migration(migrations.Migration):
    """ #2276 - the classification fields of Lab / Organization, copied into Lab/OrganizationClassificationSettings by
        classification 0233, are dropped. ClinVarKey and ClinVarKeyExcludePattern leave snpdb's state only: classification
        now owns them and their tables. """

    dependencies = [
        ('snpdb', '0290_cachedwebresource_rename_table'),
        ('classification', '0233_classification_settings'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='lab',
            name='clinvar_key',
        ),
        migrations.RemoveField(
            model_name='lab',
            name='classification_config',
        ),
        migrations.RemoveField(
            model_name='lab',
            name='consolidates_variant_classifications',
        ),
        migrations.RemoveField(
            model_name='lab',
            name='research',
        ),
        migrations.RemoveField(
            model_name='lab',
            name='upload_automatic',
        ),
        migrations.RemoveField(
            model_name='lab',
            name='upload_instructions',
        ),
        migrations.RemoveField(
            model_name='lab',
            name='upload_location',
        ),
        migrations.RemoveField(
            model_name='organization',
            name='classification_config',
        ),
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.DeleteModel(
                    name='ClinVarKeyExcludePattern',
                ),
                migrations.DeleteModel(
                    name='ClinVarKey',
                ),
            ],
            database_operations=[],
        ),
    ]
