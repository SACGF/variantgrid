# Generated by Django 4.2.10 on 2024-06-20 03:52

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0103_alter_clinvar_clinvar_review_status_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='clinvarversion',
            name='sha256_hash',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='dbnsfpgeneannotationversion',
            name='sha256_hash',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='humanproteinatlasannotationversion',
            name='sha256_hash',
            field=models.TextField(null=True),
        ),
    ]