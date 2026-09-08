from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_multiple_gnomad_sv_overlaps(apps):
    """ Only SVs that overlapped more than one gnomAD-SV record could have had the wrong one picked -
        those are the rows whose '&'-joined overlap fields hold more than one value """
    VariantAnnotation = apps.get_model("annotation", "VariantAnnotation")
    return VariantAnnotation.objects.filter(gnomad_sv_overlap_af__contains="&").exists()


class Migration(migrations.Migration):
    dependencies = [
        ("annotation", "0180_open_targets_is_lead"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["fix_gnomad_sv_overlap_chosen_record"]),
                        note="Re-pick the gnomAD-SV overlap record copied onto the gnomAD columns - "
                             "'lowest_af' compared strings, so SVs overlapping records with mixed "
                             "scientific/decimal AF notation took the most common record's values",
                        test=_has_multiple_gnomad_sv_overlaps),
    ]
