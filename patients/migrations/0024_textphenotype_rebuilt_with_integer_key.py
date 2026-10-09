""" TextPhenotype, its sentences and matches are recreated empty with an integer primary key - the sentence text was
    the key, so every sentence and match row carried a copy of it - and `match_patient_phenotypes --rebuild` splits
    and matches every description again. Every sentence awaited a rematch anyway (#2131 stamping), so there was
    nothing to keep (#2135).
    PhenotypeDescription.status was never set and DescriptionProcessingStatus had no columns, so both go. """
import django.db.models.deletion
from django.db import migrations, models

from manual.operations.manual_operations import ManualOperation


def _has_phenotype_descriptions(apps):
    PhenotypeDescription = apps.get_model("patients", "PhenotypeDescription")
    return PhenotypeDescription.objects.exists()


class Migration(migrations.Migration):

    dependencies = [
        ("ontology", "0028_backfill_phenotype_to_genes"),
        ("patients", "0023_phenotype_description_owner_links_removed"),
    ]

    operations = [
        migrations.RemoveField(model_name="phenotypedescription", name="status"),
        migrations.DeleteModel(name="DescriptionProcessingStatus"),
        migrations.DeleteModel(name="TextPhenotypeMatch"),
        migrations.DeleteModel(name="TextPhenotypeSentence"),
        migrations.DeleteModel(name="TextPhenotype"),
        migrations.CreateModel(
            name="TextPhenotype",
            fields=[
                ("id", models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("text", models.TextField(unique=True)),
                ("match_version", models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL,
                                                    to="patients.phenotypematchversion")),
            ],
        ),
        migrations.CreateModel(
            name="TextPhenotypeSentence",
            fields=[
                ("id", models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("sentence_offset", models.IntegerField()),
                ("phenotype_description", models.ForeignKey(on_delete=django.db.models.deletion.CASCADE,
                                                            to="patients.phenotypedescription")),
                ("text_phenotype", models.ForeignKey(on_delete=django.db.models.deletion.CASCADE,
                                                     to="patients.textphenotype")),
            ],
        ),
        migrations.CreateModel(
            name="TextPhenotypeMatch",
            fields=[
                ("id", models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID")),
                ("offset_start", models.IntegerField()),
                ("offset_end", models.IntegerField()),
                ("ontology_term", models.ForeignKey(on_delete=django.db.models.deletion.CASCADE,
                                                    to="ontology.ontologyterm")),
                ("text_phenotype", models.ForeignKey(on_delete=django.db.models.deletion.CASCADE,
                                                     to="patients.textphenotype")),
            ],
        ),
        ManualOperation(task_id=ManualOperation.task_id_manage(["match_patient_phenotypes", "--rebuild"]),
                        note="Split and match every phenotype description again into the integer-keyed tables (#2135)",
                        requires=["ontology-imported"],
                        test=_has_phenotype_descriptions),
    ]
