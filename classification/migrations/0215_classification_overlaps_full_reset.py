from django.db import migrations

from manual.operations.manual_operations import ManualOperation


class Migration(migrations.Migration):
    """ Builds Overlaps (with history, and the triages of old Discordance Reports) once 0211's groupings are rebuilt """

    dependencies = [
        ('classification', '0214_overlap_unique_constraints'),
    ]

    operations = [
        ManualOperation.operation_manage(["classification_overlaps", "--full_reset"],
                                         requires=[f"after:{ManualOperation.task_id_manage(['classification_groupings', '--all'])}"])
    ]
