# Generated by Django 4.1.3 on 2023-02-13 04:36

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0092_tagcolor_created_tagcolor_modified_and_more'),
        ('classification', '0093_evidence_key_data_molecular_consequence'),
    ]

    operations = [
        migrations.AlterUniqueTogether(
            name='importedalleleinfo',
            unique_together=set(),
        ),
        migrations.AddField(
            model_name='importedalleleinfo',
            name='imported_c_hgvs_md5_hash',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='importedalleleinfo',
            name='imported_g_hgvs_md5_hash',
            field=models.TextField(blank=True, null=True),
        ),
        migrations.AlterUniqueTogether(
            name='importedalleleinfo',
            unique_together={('imported_c_hgvs_md5_hash', 'imported_g_hgvs_md5_hash', 'imported_transcript', 'imported_genome_build_patch_version')},
        ),
    ]