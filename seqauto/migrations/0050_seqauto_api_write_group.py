from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _create_seqauto_api_write_group(apps, _schema_editor):
    Group = apps.get_model("auth", "Group")
    Group.objects.get_or_create(name=settings.SEQAUTO_API_WRITE_GROUP)


def _has_sequencing_runs(apps):
    SequencingRun = apps.get_model("seqauto", "SequencingRun")
    return SequencingRun.objects.exists()


class Migration(migrations.Migration):
    dependencies = [
        ("auth", "0012_alter_user_first_name_max_length"),
        ("seqauto", "0049_dragentso500combinedvariantoutput"),
    ]

    operations = [
        migrations.RunPython(_create_seqauto_api_write_group, reverse_code=migrations.RunPython.noop),
        ManualOperation.operation_other([
            f"Add the user whose API token the sequencing pipeline posts with to the "
            f"'{settings.SEQAUTO_API_WRITE_GROUP}' group (Django admin), or its writes to /seqauto/api/ get 403"
        ], test=_has_sequencing_runs),
    ]
