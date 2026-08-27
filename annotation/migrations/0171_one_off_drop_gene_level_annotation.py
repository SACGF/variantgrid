import logging

from django.db import migrations

from annotation.models.models_enums import VariantAnnotationPipelineType

# Small enough that Postgres keeps every partition on an index scan - a few hundred run ids in one
# IN list estimates high enough that the biggest partitions get sequentially scanned instead
RUN_BATCH_SIZE = 10


def _drop_gene_level_annotation(apps, _schema_editor):
    """ Gene-level annotation written through the ORM landed in the base table the version partitions
        inherit from, and every annotation queryset reads the partition directly - so grids, gene lists
        and comp-het were all blind to it. Those rows also predate the representative gene/transcript.

        Deleting the runs cascades the rows away; the scheduler then recreates a GENE_LEVEL run for
        each existing range lock and the pipeline writes them again, into the partition and with a
        gene. @see annotation.gene_level_annotation

        The cascade reaches the annotation tables by annotation_run_id, which inheritance turns into a
        delete against every version's partition, so it goes a batch of runs at a time - a whole
        deployment's worth of run ids in one statement is hundreds of GB of sequential scan. """

    AnnotationRun = apps.get_model("annotation", "AnnotationRun")
    gene_level_qs = AnnotationRun.objects.filter(pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL)
    run_ids = list(gene_level_qs.values_list("pk", flat=True))
    for i in range(0, len(run_ids), RUN_BATCH_SIZE):
        batch = run_ids[i:i + RUN_BATCH_SIZE]
        AnnotationRun.objects.filter(pk__in=batch).delete()
        logging.info("Dropped gene level annotation runs: %d/%d", i + len(batch), len(run_ids))


class Migration(migrations.Migration):
    # Each batch commits as it goes, so an interrupted delete doesn't have to start again - the
    # migration only ever deletes what's left
    atomic = False

    dependencies = [
        ('annotation', '0170_one_off_backfill_sv_conservation'),
    ]

    operations = [
        migrations.RunPython(_drop_gene_level_annotation, migrations.RunPython.noop),
    ]
