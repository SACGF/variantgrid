from django.db import migrations

from annotation.models.models_enums import VariantAnnotationPipelineType


def _drop_gene_level_annotation(apps, _schema_editor):
    """ Gene-level annotation written through the ORM landed in the base table the version partitions
        inherit from, and every annotation queryset reads the partition directly - so grids, gene lists
        and comp-het were all blind to it. Those rows also predate the representative gene/transcript.

        Deleting the runs cascades the rows away; the scheduler then recreates a GENE_LEVEL run for
        each existing range lock and the pipeline writes them again, into the partition and with a
        gene. @see annotation.gene_level_annotation """

    AnnotationRun = apps.get_model("annotation", "AnnotationRun")
    AnnotationRun.objects.filter(pipeline_type=VariantAnnotationPipelineType.GENE_LEVEL).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0170_one_off_backfill_sv_conservation'),
    ]

    operations = [
        migrations.RunPython(_drop_gene_level_annotation, migrations.RunPython.noop),
    ]
