import django.db.models.deletion
from django.db import migrations, models
from django.db.models import OuterRef, Subquery


def _backfill_bam_file_sequencing_sample(apps, _schema_editor):
    BamFile = apps.get_model("seqauto", "BamFile")
    UnalignedReads = apps.get_model("seqauto", "UnalignedReads")

    sequencing_sample = UnalignedReads.objects.filter(pk=OuterRef("unaligned_reads")).values("sequencing_sample")
    BamFile.objects.filter(unaligned_reads__isnull=False).update(
        sequencing_sample=Subquery(sequencing_sample))


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0043_experiment_experiment_name_not_blank"),
    ]

    operations = [
        migrations.AddField(
            model_name="bamfile",
            name="sequencing_sample",
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE,
                                    to="seqauto.sequencingsample"),
        ),
        migrations.RunPython(_backfill_bam_file_sequencing_sample, migrations.RunPython.noop),
        migrations.AlterField(
            model_name="bamfile",
            name="sequencing_sample",
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE,
                                    to="seqauto.sequencingsample"),
        ),
    ]
