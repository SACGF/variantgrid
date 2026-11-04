"""
Expression index on length(seq) for the SV pipeline routing subquery (alt__in Sequence where Length(seq) >= N, #2635).
Built concurrently so VCF import can keep inserting sequences, then ANALYZE: the index's statistics are what let the
planner estimate the filter and keep the subplan hashed.
"""
from django.contrib.postgres.operations import AddIndexConcurrently
from django.db import migrations, models
from django.db.models.functions.text import Length


class Migration(migrations.Migration):
    atomic = False

    dependencies = [
        ("snpdb", "0006_one_off_dedupe_cohort_genotypes"),
    ]

    operations = [
        AddIndexConcurrently(
            model_name="sequence",
            index=models.Index(Length("seq"), name="snpdb_sequence_seq_length_idx"),
        ),
        migrations.RunSQL("ANALYZE snpdb_sequence", migrations.RunSQL.noop),
    ]
