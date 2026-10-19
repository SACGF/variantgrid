import django.db.models.deletion
import django_extensions.db.fields
from django.db import migrations, models


def _copy_settings_from_lab_and_organization(apps, _schema_editor):
    Organization = apps.get_model("snpdb", "Organization")
    Lab = apps.get_model("snpdb", "Lab")
    OrganizationClassificationSettings = apps.get_model("classification", "OrganizationClassificationSettings")
    LabClassificationSettings = apps.get_model("classification", "LabClassificationSettings")

    OrganizationClassificationSettings.objects.bulk_create([
        OrganizationClassificationSettings(organization_id=org_id, classification_config=classification_config)
        for org_id, classification_config in Organization.objects.values_list("pk", "classification_config")
    ])
    LabClassificationSettings.objects.bulk_create([
        LabClassificationSettings(lab_id=lab.pk,
                                  clinvar_key_id=lab.clinvar_key_id,
                                  classification_config=lab.classification_config,
                                  research=lab.research,
                                  consolidates_variant_classifications=lab.consolidates_variant_classifications,
                                  upload_location=lab.upload_location,
                                  upload_automatic=lab.upload_automatic,
                                  upload_instructions=lab.upload_instructions)
        for lab in Lab.objects.all()
    ])


def _copy_settings_back_to_lab_and_organization(apps, _schema_editor):
    Organization = apps.get_model("snpdb", "Organization")
    Lab = apps.get_model("snpdb", "Lab")
    OrganizationClassificationSettings = apps.get_model("classification", "OrganizationClassificationSettings")
    LabClassificationSettings = apps.get_model("classification", "LabClassificationSettings")

    for org_settings in OrganizationClassificationSettings.objects.all():
        Organization.objects.filter(pk=org_settings.organization_id).update(
            classification_config=org_settings.classification_config)
    for lab_settings in LabClassificationSettings.objects.all():
        Lab.objects.filter(pk=lab_settings.lab_id).update(
            clinvar_key_id=lab_settings.clinvar_key_id,
            classification_config=lab_settings.classification_config,
            research=lab_settings.research,
            consolidates_variant_classifications=lab_settings.consolidates_variant_classifications,
            upload_location=lab_settings.upload_location,
            upload_automatic=lab_settings.upload_automatic,
            upload_instructions=lab_settings.upload_instructions)


class Migration(migrations.Migration):
    """ #2276 - ClinVarKey and ClinVarKeyExcludePattern move from snpdb in state only, keeping their tables, and the
        classification fields of Lab / Organization are copied into the new Lab/OrganizationClassificationSettings.
        snpdb 0291 then drops those fields. """

    dependencies = [
        ('classification', '0232_overlapdiscordancenotification_lab'),
        ('snpdb', '0290_cachedwebresource_rename_table'),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.CreateModel(
                    name='ClinVarKey',
                    fields=[
                        ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                        ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                        ('id', models.TextField(primary_key=True, serialize=False)),
                        ('api_key', models.TextField(blank=True, null=True)),
                        ('org_id', models.TextField(blank=True, default='')),
                        ('name', models.TextField(blank=True, default='')),
                        ('last_full_run', models.DateTimeField(null=True)),
                        ('default_affected_status', models.TextField(blank=True, choices=[('yes', 'Yes'), ('no', 'No'), ('unknown', 'Unknown'), ('not provided', 'Not Provided'), ('not applicable', 'Not Applicable')], null=True)),
                        ('inject_acmg_description', models.BooleanField(blank=True, default=False)),
                        ('include_interpretation_summary', models.BooleanField(blank=True, default=True)),
                        ('assertion_method_lookup', models.JSONField(default=dict)),
                        ('citations_mode', models.TextField(choices=[('all', 'All'), ('interpret', 'Interpretation Summary Only')], default='all')),
                    ],
                    options={
                        'verbose_name': 'ClinVar key',
                        'db_table': 'snpdb_clinvarkey',
                    },
                ),
                migrations.AlterField(
                    model_name='clinvarallele',
                    name='clinvar_key',
                    field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='classification.clinvarkey'),
                ),
                migrations.AlterField(
                    model_name='clinvarexportbatch',
                    name='clinvar_key',
                    field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='classification.clinvarkey'),
                ),
                migrations.CreateModel(
                    name='ClinVarKeyExcludePattern',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                        ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                        ('evidence_key', models.TextField()),
                        ('pattern', models.TextField()),
                        ('name', models.TextField(blank=True)),
                        ('case_insensitive', models.BooleanField(blank=True, default=True)),
                        ('mode', models.TextField(choices=[('X', 'Exclude If Match'), ('N', 'Exclude If No Match')], default='X')),
                        ('clinvar_key', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='classification.clinvarkey')),
                    ],
                    options={
                        'db_table': 'snpdb_clinvarkeyexcludepattern',
                    },
                ),
            ],
            database_operations=[],
        ),
        migrations.CreateModel(
            name='OrganizationClassificationSettings',
            fields=[
                ('organization', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, related_name='classification_settings', serialize=False, to='snpdb.organization')),
                ('classification_config', models.JSONField(blank=True, null=True)),
            ],
            options={
                'verbose_name_plural': 'organization classification settings',
                'base_manager_name': 'objects',
            },
        ),
        migrations.CreateModel(
            name='LabClassificationSettings',
            fields=[
                ('lab', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, primary_key=True, related_name='classification_settings', serialize=False, to='snpdb.lab')),
                ('classification_config', models.JSONField(blank=True, null=True)),
                ('research', models.BooleanField(blank=True, default=False)),
                ('consolidates_variant_classifications', models.BooleanField(default=False)),
                ('upload_location', models.TextField(blank=True, null=True)),
                ('upload_automatic', models.BooleanField(blank=True, default=False)),
                ('upload_instructions', models.TextField(blank=True, default='')),
                ('clinvar_key', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='lab_classification_settings', to='classification.clinvarkey')),
            ],
            options={
                'verbose_name_plural': 'lab classification settings',
                'base_manager_name': 'objects',
            },
        ),
        migrations.RunPython(_copy_settings_from_lab_and_organization, _copy_settings_back_to_lab_and_organization),
    ]
