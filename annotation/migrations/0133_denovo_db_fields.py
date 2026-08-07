from django.db import migrations, models


def _new_variant_annotation_fields():
    return [
        ('denovo_db_studies', models.TextField(blank=True, null=True)),
        ('denovo_db_pubmed_ids', models.TextField(blank=True, null=True)),
        ('denovo_db_primary_phenotypes', models.TextField(blank=True, null=True)),
        ('denovo_db_case_count', models.IntegerField(blank=True, null=True)),
        ('denovo_db_control_count', models.IntegerField(blank=True, null=True)),
    ]


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0132_alter_variantannotationversion_active'),
    ]

    operations = [
        migrations.AddField(
            model_name='variantannotationversion',
            name='denovo_db',
            field=models.TextField(blank=True, null=True),
        ),
        *[
            migrations.AddField(model_name='variantannotation', name=name, field=field)
            for name, field in _new_variant_annotation_fields()
        ],
    ]
