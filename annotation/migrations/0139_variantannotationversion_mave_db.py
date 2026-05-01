import datetime

from django.db import migrations, models
from django.utils import timezone


# MaveDB plugin was added with the file MaveDB_variants_2023-11-29.tsv.gz; any
# VariantAnnotationVersion created after that date would have included MaveDB.
MAVE_DB_INITIAL_VERSION = "2023-11-29"


def _backfill_mave_db(apps, _schema_editor):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    cutoff = timezone.make_aware(datetime.datetime(2023, 11, 29))
    VariantAnnotationVersion.objects.filter(
        annotation_date__gt=cutoff,
        genome_build_id="GRCh38",
    ).update(mave_db=MAVE_DB_INITIAL_VERSION)


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0138_cohort_genotype_annotation_stats_backfill"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="mave_db",
            field=models.TextField(blank=True, null=True),
        ),
        migrations.RunPython(_backfill_mave_db, migrations.RunPython.noop),
    ]
