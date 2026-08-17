from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('upload', '0032_one_off_backfill_uploaded_vcf_pipeline_max_variant'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='UploadedFile',
            new_name='FileUpload',
        ),
        migrations.RenameField(
            model_name='fileupload',
            old_name='uploaded_file',
            new_name='file_field',
        ),
        migrations.RenameField(
            model_name='uploadpipeline',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedvcf',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedbed',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedclassificationimport',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedclinvarversion',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedgenecoverage',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedgenelist',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedliftover',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedmanualvariantentrycollection',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedpatientrecords',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedpedfile',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedvarianttags',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
        migrations.RenameField(
            model_name='uploadedwikicollection',
            old_name='uploaded_file',
            new_name='file_upload',
        ),
    ]
