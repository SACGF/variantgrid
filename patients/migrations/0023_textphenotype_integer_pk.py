""" TextPhenotype gets an integer primary key - the sentence text was the key, so every sentence and match row
    carried a copy of it (#2135). Rows are numbered in text order and the two FKs repointed by joining on text. """
import django.db.models.deletion
from django.db import migrations, models

NUMBER_TEXT_PHENOTYPES_SQL = """
UPDATE patients_textphenotype t SET id = n.row_number
FROM (SELECT text, row_number() OVER (ORDER BY text) FROM patients_textphenotype) n
WHERE n.text = t.text
"""

# Django adds the identity without moving it past the numbered rows
SET_IDENTITY_SQL = """
SELECT setval(pg_get_serial_sequence('patients_textphenotype', 'id'), coalesce(max(id), 1), max(id) IS NOT NULL)
FROM patients_textphenotype
"""

REPOINT_SQL = [
    f"""UPDATE {table} r SET text_phenotype_new = t.id
        FROM patients_textphenotype t WHERE t.text = r.text_phenotype_id"""
    for table in ["patients_textphenotypesentence", "patients_textphenotypematch"]
] + [
    # Check the deferred FKs now, or the ALTER TABLEs below fail on pending trigger events
    "SET CONSTRAINTS ALL IMMEDIATE",
]


class Migration(migrations.Migration):

    dependencies = [
        ("patients", "0022_remove_phenotype_status_and_processed"),
    ]

    operations = [
        migrations.AddField(
            model_name="textphenotype",
            name="id",
            field=models.IntegerField(null=True),
        ),
        migrations.RunSQL(NUMBER_TEXT_PHENOTYPES_SQL),
        migrations.AddField(
            model_name="textphenotypesentence",
            name="text_phenotype_new",
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name="textphenotypematch",
            name="text_phenotype_new",
            field=models.IntegerField(null=True),
        ),
        migrations.RunSQL(REPOINT_SQL),
        migrations.RemoveField(model_name="textphenotypesentence", name="text_phenotype"),
        migrations.RemoveField(model_name="textphenotypematch", name="text_phenotype"),
        migrations.AlterField(
            model_name="textphenotype",
            name="text",
            field=models.TextField(unique=True),
        ),
        migrations.AlterField(
            model_name="textphenotype",
            name="id",
            field=models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name="ID"),
        ),
        migrations.RunSQL(SET_IDENTITY_SQL),
        migrations.AlterField(
            model_name="textphenotypesentence",
            name="text_phenotype_new",
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name="+",
                                    to="patients.textphenotype"),
        ),
        migrations.AlterField(
            model_name="textphenotypematch",
            name="text_phenotype_new",
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name="+",
                                    to="patients.textphenotype"),
        ),
        migrations.RenameField(model_name="textphenotypesentence", old_name="text_phenotype_new",
                               new_name="text_phenotype"),
        migrations.RenameField(model_name="textphenotypematch", old_name="text_phenotype_new",
                               new_name="text_phenotype"),
        migrations.AlterField(
            model_name="textphenotypesentence",
            name="text_phenotype",
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to="patients.textphenotype"),
        ),
        migrations.AlterField(
            model_name="textphenotypematch",
            name="text_phenotype",
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to="patients.textphenotype"),
        ),
    ]
