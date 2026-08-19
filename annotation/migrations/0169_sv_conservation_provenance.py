from django.db import migrations, models

# The --custom short_names an SV VEP command carried while VEP itself scored conservation. Frozen here
# rather than read from the VEPCustom enum: this migration classifies the runs that exist right now.
CONSERVATION_SHORT_NAMES = [
    "phastCons100way", "phastCons30way", "phastCons46way",
    "phyloP100way", "phyloP30way", "phyloP46way",
]

STRUCTURAL_VARIANT = "C"  # VariantAnnotationPipelineType.STRUCTURAL_VARIANT


def _pybigwig_scored_runs(AnnotationRun):
    """ SV runs whose VEP command carried no conservation --custom args, ie the ones where our pyBigWig
        stage was responsible for those columns (#1657). The command is the only record of that until
        sv_conservation_pybigwig starts being set at annotate time. """
    qs = AnnotationRun.objects.filter(pipeline_type=STRUCTURAL_VARIANT).exclude(pipeline_command=None)
    for short_name in CONSERVATION_SHORT_NAMES:
        qs = qs.exclude(pipeline_command__contains=f"short_name={short_name}")
    return qs


def _record_sv_conservation_provenance(apps, _schema_editor):
    AnnotationRun = apps.get_model("annotation", "AnnotationRun")
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")

    run_qs = _pybigwig_scored_runs(AnnotationRun)
    version_ids = set(run_qs.values_list("annotation_range_lock__version_id", flat=True))
    run_qs.update(sv_conservation_pybigwig=True)

    # Those versions carry the values the fixed window and the failed sidecars have to be re-scored into -
    # but only where conservation is annotated at all. A build with no conservation bigWigs (T2T) has
    # pyBigWig-scored runs too, and nothing for the backfill to put in them.
    needs_backfill = [version_id for version_id in version_ids
                      if VariantAnnotation.objects.filter(version_id=version_id,
                                                          phastcons_100_way_vertebrate__isnull=False).exists()]
    VariantAnnotationVersion.objects.filter(pk__in=needs_backfill).update(backfilled_sv_conservation=False)


def _clear_sv_conservation_provenance(apps, _schema_editor):
    apps.get_model("annotation", "AnnotationRun").objects.update(sv_conservation_pybigwig=False)
    apps.get_model("annotation", "VariantAnnotationVersion").objects.update(backfilled_sv_conservation=True)


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0168_annotation_pipeline_version_lifecycle'),
    ]

    operations = [
        migrations.AddField(
            model_name='annotationrun',
            name='sv_conservation_pybigwig',
            field=models.BooleanField(default=False),
        ),
        migrations.AddField(
            model_name='variantannotationversion',
            name='backfilled_sv_conservation',
            field=models.BooleanField(default=True),
        ),
        migrations.RunPython(_record_sv_conservation_provenance, _clear_sv_conservation_provenance),
    ]
