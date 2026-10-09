""" The phenotype text models move from annotation to patients (#2135): the tables are renamed and patients' state
    takes the models exactly as annotation had them; annotation 0193 then drops them from its state.
    Index, constraint and sequence names keep their annotation_ prefix. """
import django.db.models.deletion
import django.utils.timezone
import model_utils.fields
from django.conf import settings
from django.db import migrations, models

PHENOTYPE_TABLES = [
    "descriptionprocessingstatus",
    "phenotypedescription",
    "phenotypematchversion",
    "textphenotype",
    "textphenotypesentence",
    "textphenotypematch",
    "patienttextphenotype",
    "cohorttextphenotype",
]


def _rename_tables(old_app: str, new_app: str) -> list[str]:
    return [f"ALTER TABLE {old_app}_{table} RENAME TO {new_app}_{table}" for table in PHENOTYPE_TABLES]


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0192_one_off_match_patient_phenotypes_stale_fuzzy"),
        ("ontology", "0028_backfill_phenotype_to_genes"),
        ("patients", "0019_specimen_measure_pathology_only"),
        ("snpdb", "0281_genomic_intervals_collection_error_message"),
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.SeparateDatabaseAndState(
            database_operations=[
                migrations.RunSQL(_rename_tables("annotation", "patients"),
                                  reverse_sql=_rename_tables("patients", "annotation")),
            ],
            state_operations=[
                migrations.CreateModel(
                    name='DescriptionProcessingStatus',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                    ],
                ),
                migrations.CreateModel(
                    name='PhenotypeDescription',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('original_text', models.TextField()),
                        ('status', models.CharField(choices=[('C', 'Created'), ('T', 'Tokenized'), ('E', 'Error'), ('S', 'Success')], max_length=1)),
                    ],
                ),
                migrations.CreateModel(
                    name='PhenotypeMatchVersion',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name='created')),
                        ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name='modified')),
                        ('matcher_version', models.IntegerField()),
                        ('ontology_version', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, to='ontology.ontologyversion')),
                    ],
                    options={
                        'constraints': [models.UniqueConstraint(fields=('matcher_version', 'ontology_version'), name='phenotype_match_version_unique', nulls_distinct=False)],
                    },
                ),
                migrations.CreateModel(
                    name='TextPhenotype',
                    fields=[
                        ('text', models.TextField(primary_key=True, serialize=False)),
                        ('processed', models.BooleanField(default=False)),
                        ('match_version', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, to='patients.phenotypematchversion')),
                    ],
                ),
                migrations.CreateModel(
                    name='TextPhenotypeSentence',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('sentence_offset', models.IntegerField()),
                        ('phenotype_description', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.phenotypedescription')),
                        ('text_phenotype', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.textphenotype')),
                    ],
                ),
                migrations.CreateModel(
                    name='TextPhenotypeMatch',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('offset_start', models.IntegerField()),
                        ('offset_end', models.IntegerField()),
                        ('text_phenotype', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='patients.textphenotype')),
                        ('ontology_term', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ontology.ontologyterm')),
                    ],
                ),
                migrations.CreateModel(
                    name='PatientTextPhenotype',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('approved_by', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to=settings.AUTH_USER_MODEL)),
                        ('patient', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='patient_text_phenotype', to='patients.patient')),
                        ('phenotype_description', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='patients.phenotypedescription')),
                    ],
                ),
                migrations.CreateModel(
                    name='CohortTextPhenotype',
                    fields=[
                        ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                        ('approved_by', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, to=settings.AUTH_USER_MODEL)),
                        ('cohort', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, related_name='cohort_text_phenotype', to='snpdb.cohort')),
                        ('phenotype_description', models.OneToOneField(on_delete=django.db.models.deletion.CASCADE, to='patients.phenotypedescription')),
                    ],
                ),
            ],
        ),
    ]
