import django.db.models.deletion
from django.db import migrations, models


def _copy_alignment_file_to_m2m_and_sequencing_sample(apps, schema_editor):
    AlignmentFile = apps.get_model("seqauto", "AlignmentFile")
    sequencing_sample_by_alignment_file = dict(AlignmentFile.objects.values_list("pk", "sequencing_sample_id"))

    for model_name in ["SingleSampleVCF", "QC"]:
        klass = apps.get_model("seqauto", model_name)
        through = klass.alignment_files.through
        through_fk = f"{model_name.lower()}_id"
        records = []
        through_records = []
        for record in klass.objects.all().only("pk", "alignment_file_id"):
            record.sequencing_sample_id = sequencing_sample_by_alignment_file[record.alignment_file_id]
            records.append(record)
            through_records.append(through(**{through_fk: record.pk, "alignmentfile_id": record.alignment_file_id}))
        klass.objects.bulk_update(records, ["sequencing_sample"], batch_size=2000)
        through.objects.bulk_create(through_records, batch_size=2000)
    # Deferred FK checks would otherwise block the ALTER TABLEs below in the same transaction
    schema_editor.execute("SET CONSTRAINTS ALL IMMEDIATE")


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0052_alignment_file"),
    ]

    operations = [
        migrations.AddField(
            model_name="singlesamplevcf",
            name="sequencing_sample",
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE,
                                    to="seqauto.sequencingsample"),
        ),
        migrations.AddField(
            model_name="singlesamplevcf",
            name="alignment_files",
            field=models.ManyToManyField(to="seqauto.alignmentfile"),
        ),
        migrations.AddField(
            model_name="qc",
            name="sequencing_sample",
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE,
                                    to="seqauto.sequencingsample"),
        ),
        migrations.AddField(
            model_name="qc",
            name="alignment_files",
            field=models.ManyToManyField(to="seqauto.alignmentfile"),
        ),
        migrations.RunPython(_copy_alignment_file_to_m2m_and_sequencing_sample, migrations.RunPython.noop),
        migrations.AlterField(
            model_name="singlesamplevcf",
            name="sequencing_sample",
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to="seqauto.sequencingsample"),
        ),
        migrations.AlterField(
            model_name="qc",
            name="sequencing_sample",
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to="seqauto.sequencingsample"),
        ),
        migrations.RemoveField(
            model_name="singlesamplevcf",
            name="alignment_file",
        ),
        migrations.RemoveField(
            model_name="qc",
            name="alignment_file",
        ),
    ]
