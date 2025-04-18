# Generated by Django 4.0.3 on 2022-04-27 05:52

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('upload', '0013_alter_uploadedfile_file_type_uploadedwikicollection'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='uploadedfile',
            name='visible',
        ),
        migrations.RemoveField(
            model_name='uploadsettings',
            name='show_all',
        ),
        migrations.CreateModel(
            name='UploadSettingsFileType',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('file_type', models.CharField(choices=[('B', 'BED'), ('L', 'Clinvar'), ('C', 'CuffDiff'), ('G', 'Gene List'), ('O', 'Gene Coverage'), ('I', 'Liftover'), ('P', 'Pedigree'), ('R', 'Patient Records'), ('S', 'Variant Classifications'), ('A', 'Variant Tags'), ('V', 'VCF'), ('Y', 'VCF - Insert variants only (no samples etc)'), ('w', 'Gene Wiki records'), ('W', 'Variant Wiki records')], max_length=1, null=True)),
                ('upload_settings', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='upload.uploadsettings')),
            ],
            options={
                'unique_together': {('upload_settings', 'file_type')},
            },
        ),
    ]
