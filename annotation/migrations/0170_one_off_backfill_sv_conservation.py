from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_sv_conservation_to_backfill(apps):
    """ 0169 flips this False on every annotation version whose SV conservation our pyBigWig stage
        wrote - those are the values the corrected window has to be re-scored into. """
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(backfilled_sv_conservation=False).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0169_sv_conservation_provenance"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["backfill_sv_conservation"]),
                        note="Recompute the SV conservation (phastCons/phyloP max) columns - the pyBigWig "
                             "scoring window started one base early, and a run whose sidecar stage failed "
                             "left the columns null (#1657)",
                        test=_has_sv_conservation_to_backfill),
    ]
