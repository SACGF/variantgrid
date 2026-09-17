from django.db import migrations


def _rename_gene_level_api_uploads(apps, _schema_editor):
    """ Both classification import pipelines used to be called 'Variants from API', so the upload
        listing showed two rows with the same name. @see classification.classification_import """
    FileUpload = apps.get_model("upload", "FileUpload")
    FileUpload.objects.filter(file_type='y', name='Variants from API').update(name='Gene-level Variants from API')


class Migration(migrations.Migration):
    dependencies = [
        ('upload', '0044_alter_fileupload_file_type_and_more'),
    ]

    operations = [
        migrations.RunPython(_rename_gene_level_api_uploads, migrations.RunPython.noop),
    ]
