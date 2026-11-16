from django.db import migrations, models
from django.db.models import IntegerField, Max, Min
from django.db.models.fields.json import KeyTextTransform, KeyTransform
from django.db.models.functions import Cast

BATCH_SIZE = 50_000


def _backfill_transcript_version_spans(apps, _schema_editor):
    """ start/end/strand out of each row's cdot data for its build, in SQL - the rows are never loaded """
    TranscriptVersion = apps.get_model("genes", "TranscriptVersion")
    for genome_build_id in TranscriptVersion.objects.values_list("genome_build_id", flat=True).distinct().order_by():
        build_data = KeyTransform(genome_build_id, KeyTransform("genome_builds", "data"))
        exons = KeyTransform("exons", build_data)
        build_qs = TranscriptVersion.objects.filter(genome_build_id=genome_build_id)
        pk_range = build_qs.aggregate(min_pk=Min("pk"), max_pk=Max("pk"))
        for batch_start in range(pk_range["min_pk"], pk_range["max_pk"] + 1, BATCH_SIZE):
            batch_qs = build_qs.filter(pk__gte=batch_start, pk__lt=batch_start + BATCH_SIZE)
            batch_qs.update(start=Cast(KeyTextTransform(0, KeyTransform(0, exons)), IntegerField()),
                            end=Cast(KeyTextTransform(1, KeyTransform(-1, exons)), IntegerField()),
                            strand=KeyTextTransform("strand", build_data))


class Migration(migrations.Migration):

    dependencies = [
        ("genes", "0005_cdot_data_version_transcript_stats"),
    ]

    operations = [
        migrations.AddField(
            model_name="transcriptversion",
            name="start",
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name="transcriptversion",
            name="end",
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name="transcriptversion",
            name="strand",
            field=models.CharField(max_length=1, null=True),
        ),
        migrations.RunPython(_backfill_transcript_version_spans, migrations.RunPython.noop),
        migrations.AddIndex(
            model_name="transcriptversion",
            index=models.Index(fields=["contig", "start"], include=["end", "gene_version", "genome_build", "id"],
                               name="genes_tv_contig_span_idx"),
        ),
    ]
