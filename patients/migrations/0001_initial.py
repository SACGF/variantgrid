# Generated by Django 3.1 on 2020-09-29 05:30

import django.db.models.deletion
import django_extensions.db.fields
from django.conf import settings
from django.db import migrations, models

import annotation.models.has_phenotype_description_mixin
import library.django_utils.django_file_system_storage
import library.django_utils.guardian_permissions_mixin


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='ExternalModelManager',
            fields=[
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('name', models.TextField(primary_key=True, serialize=False)),
                ('details', models.TextField(blank=True)),
                ('can_modify', models.BooleanField(default=False)),
                ('modifications_sent_to_external_system', models.BooleanField(default=False)),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='ExternalPK',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('code', models.TextField()),
                ('external_type', models.TextField()),
                ('external_manager', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.externalmodelmanager')),
            ],
            options={
                'unique_together': {('code', 'external_type', 'external_manager')},
            },
        ),
        migrations.CreateModel(
            name='FakeData',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
            ],
            options={
                'get_latest_by': 'modified',
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Patient',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('family_code', models.TextField(blank=True, null=True)),
                ('first_name', models.TextField(blank=True, null=True)),
                ('last_name', models.TextField(null=True)),
                ('date_of_birth', models.DateField(blank=True, null=True)),
                ('date_of_death', models.DateField(blank=True, null=True)),
                ('sex', models.CharField(choices=[('U', 'unknown'), ('M', 'male'), ('F', 'female')], default='U', max_length=1)),
                ('phenotype', models.TextField(blank=True, null=True)),
                ('affected', models.BooleanField(null=True)),
                ('consanguineous', models.BooleanField(null=True)),
                ('medicare', models.TextField(blank=True, null=True)),
                ('billing_details', models.TextField(blank=True, null=True)),
                ('street_address', models.TextField(blank=True, null=True)),
                ('suburb', models.TextField(blank=True, null=True)),
                ('postcode', models.TextField(blank=True, null=True)),
                ('state', models.TextField(blank=True, null=True)),
                ('telephone', models.IntegerField(blank=True, null=True)),
                ('_deceased', models.BooleanField(blank=True, null=True)),
                ('external_pk', models.OneToOneField(null=True, on_delete=django.db.models.deletion.CASCADE, to='patients.externalpk')),
                ('fake_data', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='patients.fakedata')),
            ],
            options={
                'abstract': False,
            },
            bases=(library.django_utils.guardian_permissions_mixin.GuardianPermissionsMixin, annotation.models.has_phenotype_description_mixin.HasPhenotypeDescriptionMixin, models.Model),
        ),
        migrations.CreateModel(
            name='PatientImport',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='Tissue',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.TextField()),
                ('description', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='Specimen',
            fields=[
                ('reference_id', models.TextField(primary_key=True, serialize=False)),
                ('description', models.TextField(blank=True, null=True)),
                ('collected_by', models.TextField(blank=True, null=True)),
                ('collection_date', models.DateTimeField(blank=True, null=True)),
                ('received_date', models.DateTimeField(blank=True, null=True)),
                ('mutation_type', models.CharField(blank=True, choices=[('G', 'Germline'), ('S', 'Somatic')], default='G', max_length=1, null=True)),
                ('nucleic_acid_source', models.CharField(blank=True, choices=[('D', 'DNA'), ('R', 'RNA')], default='D', max_length=1, null=True)),
                ('_age_at_collection_date', models.IntegerField(blank=True, null=True)),
                ('patient', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.patient')),
                ('tissue', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='patients.tissue')),
            ],
        ),
        migrations.CreateModel(
            name='PatientRecords',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('patient_import', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='patients.patientimport')),
            ],
        ),
        migrations.CreateModel(
            name='PatientRecord',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('record_id', models.IntegerField()),
                ('valid', models.BooleanField(default=False)),
                ('validation_message', models.TextField(blank=True, null=True)),
                ('matched_sample_id', models.IntegerField(null=True)),
                ('sample_id', models.IntegerField(null=True)),
                ('sample_name', models.TextField(null=True)),
                ('patient_family_code', models.TextField(null=True)),
                ('patient_first_name', models.TextField(null=True)),
                ('patient_last_name', models.TextField()),
                ('date_of_birth', models.DateField(null=True)),
                ('date_of_death', models.DateField(null=True)),
                ('sex', models.CharField(choices=[('U', 'unknown'), ('M', 'male'), ('F', 'female')], max_length=1, null=True)),
                ('affected', models.BooleanField(blank=True, null=True)),
                ('consanguineous', models.BooleanField(blank=True, null=True)),
                ('_deceased', models.BooleanField(blank=True, null=True)),
                ('patient_phenotype', models.TextField(null=True)),
                ('specimen_reference_id', models.TextField(null=True)),
                ('specimen_description', models.TextField(null=True)),
                ('specimen_collected_by', models.TextField(null=True)),
                ('specimen_collection_date', models.TextField(null=True)),
                ('specimen_received_date', models.TextField(null=True)),
                ('specimen_mutation_type', models.CharField(choices=[('G', 'Germline'), ('S', 'Somatic')], max_length=1, null=True)),
                ('specimen_nucleic_acid_source', models.CharField(choices=[('D', 'DNA'), ('R', 'RNA')], max_length=1, null=True)),
                ('specimen_age_at_collection_date', models.IntegerField(blank=True, null=True)),
                ('created_patient', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='created_patient', to='patients.patient')),
                ('created_specimen', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='created_specimen', to='patients.specimen')),
                ('matched_patient', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='matched_patient', to='patients.patient')),
                ('matched_specimen', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='matched_specimen', to='patients.specimen')),
                ('patient_records', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.patientrecords')),
            ],
        ),
        migrations.CreateModel(
            name='PatientPopulation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('population', models.CharField(choices=[('AFR', 'African/African American'), ('ASJ', 'Ashkenazi Jewish'), ('AM', 'Australo Melanesian'), ('CA', 'Central Asian'), ('EAS', 'East Asian'), ('FIN', 'Finnish'), ('AMR', 'Latino / Mixed Amerindian'), ('MEA', 'Middle East / North African'), ('NFE', 'Non-Finnish European'), ('PO', 'Polynesian'), ('SAS', 'South Asian'), ('SEA', 'South East Asian'), ('OTH', 'Other')], max_length=3)),
                ('patient', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.patient')),
            ],
        ),
        migrations.CreateModel(
            name='PatientModification',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('date', models.DateTimeField(auto_now_add=True)),
                ('description', models.TextField(null=True)),
                ('origin', models.CharField(choices=[('C', 'Uploaded CSV'), ('I', 'Internal change'), ('G', 'Manual change by user'), ('E', 'External Database')], max_length=1)),
                ('patient', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.patient')),
                ('patient_import', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='patients.patientimport')),
                ('user', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='PatientComment',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('comment', models.TextField()),
                ('patient', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.patient')),
                ('patient_modification', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.patientmodification')),
            ],
        ),
        migrations.CreateModel(
            name='PatientAttachment',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file', models.FileField(storage=library.django_utils.django_file_system_storage.PrivateUploadStorage(), upload_to='patient_attachments')),
                ('file_type', models.CharField(choices=[('I', 'Image'), ('P', 'Image'), ('W', 'Word Doc'), ('T', 'Text'), ('S', 'Spreadsheet'), ('O', 'Image')], max_length=1)),
                ('thumbnail_path', models.TextField(null=True)),
                ('patient', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.patient')),
            ],
        ),
        migrations.CreateModel(
            name='FollowLeadScientist',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('follow', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='followed', to=settings.AUTH_USER_MODEL)),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='following', to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='Clinician',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('email', models.TextField(blank=True, null=True)),
                ('title', models.CharField(blank=True, choices=[('A', 'Associate Professor'), ('D', 'Dr'), ('R', 'Mr'), ('S', 'Mrs'), ('M', 'Ms'), ('P', 'Professor')], max_length=1, null=True)),
                ('first_name', models.TextField(blank=True, null=True)),
                ('last_name', models.TextField(null=True)),
                ('specialty', models.TextField(blank=True, null=True)),
                ('phone', models.TextField(blank=True, null=True)),
                ('user', models.OneToOneField(null=True, on_delete=django.db.models.deletion.SET_NULL, to=settings.AUTH_USER_MODEL)),
            ],
        ),
    ]
