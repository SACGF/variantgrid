from django.db import migrations


def _cdot_data_version_from_transcripts(apps, _schema_editor):
    """ Before #2029 every import rewrote all of a build/consortium's TranscriptVersions, so the most recently
        inserted one carries the cdot version that's installed """
    CdotDataVersion = apps.get_model("genes", "CdotDataVersion")
    GeneAnnotationImport = apps.get_model("genes", "GeneAnnotationImport")
    TranscriptVersion = apps.get_model("genes", "TranscriptVersion")

    build_and_consortia = GeneAnnotationImport.objects.values_list("genome_build_id", "annotation_consortium").distinct()
    for genome_build_id, annotation_consortium in build_and_consortia:
        tv_qs = TranscriptVersion.objects.filter(genome_build_id=genome_build_id,
                                                 transcript__annotation_consortium=annotation_consortium)
        cdot_version = tv_qs.order_by("-pk").values_list("data__cdot", flat=True).first()
        if cdot_version:
            CdotDataVersion.objects.update_or_create(genome_build_id=genome_build_id,
                                                     annotation_consortium=annotation_consortium,
                                                     defaults={"cdot_version": cdot_version})


class Migration(migrations.Migration):
    dependencies = [
        ("genes", "0097_cdot_data_version"),
    ]

    operations = [
        migrations.RunPython(_cdot_data_version_from_transcripts, migrations.RunPython.noop),
    ]
