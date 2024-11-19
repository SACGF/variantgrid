from django.db import migrations

from manual.operations.manual_operations import ManualOperation


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0114_classification_withdraw_reason'),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["remove_invisible_characters"]))
    ]
