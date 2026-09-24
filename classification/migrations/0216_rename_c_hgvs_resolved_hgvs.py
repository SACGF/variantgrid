from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("classification", "0215_classification_overlaps_full_reset"),
    ]

    operations = [
        migrations.RenameField(
            model_name="resolvedvariantinfo",
            old_name="c_hgvs",
            new_name="resolved_hgvs",
        ),
        migrations.RenameField(
            model_name="resolvedvariantinfo",
            old_name="c_hgvs_compat",
            new_name="resolved_hgvs_compat",
        ),
        migrations.RenameField(
            model_name="resolvedvariantinfo",
            old_name="c_hgvs_converter_version",
            new_name="hgvs_converter_version",
        ),
        migrations.RenameField(
            model_name="resolvedvariantinfo",
            old_name="c_hgvs_converter_data_version",
            new_name="hgvs_converter_data_version",
        ),
        migrations.RenameField(
            model_name="importedalleleinfovalidation",
            old_name="c_hgvs_37",
            new_name="resolved_hgvs_37",
        ),
        migrations.RenameField(
            model_name="importedalleleinfovalidation",
            old_name="c_hgvs_38",
            new_name="resolved_hgvs_38",
        ),
    ]
