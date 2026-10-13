from django.db import migrations

# QC lost its single alignment_file FK (0053_alignment_files_m2m: alignment_files M2M plus a direct sequencing_sample),
# so QCType.qc_object_path rows that 0053_qctype_alignment_file_path pointed at alignment_file__ need the new routes.
# Longest prefix first: each path matches exactly one.
FORWARD_PREFIXES = [
    ("alignment_file__unaligned_reads__sequencing_sample__", "sequencing_sample__"),
    ("alignment_file__unaligned_reads__", "sequencing_sample__unalignedreads__"),
    ("alignment_file__", "alignment_files__"),
]


def _rewrite(apps, prefixes):
    QCType = apps.get_model("seqauto", "QCType")
    for qc_type in QCType.objects.all():
        for old, new in prefixes:
            if qc_type.qc_object_path.startswith(old):
                qc_type.qc_object_path = new + qc_type.qc_object_path[len(old):]
                qc_type.save(update_fields=["qc_object_path"])
                break


def _forward(apps, _schema_editor):
    _rewrite(apps, FORWARD_PREFIXES)


def _reverse(apps, _schema_editor):
    _rewrite(apps, [(new, old) for old, new in FORWARD_PREFIXES])


class Migration(migrations.Migration):

    dependencies = [
        ("seqauto", "0053_alignment_files_m2m"),
        ("seqauto", "0053_qctype_alignment_file_path"),
    ]

    operations = [
        migrations.RunPython(_forward, _reverse),
    ]
