from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _user_awards_enabled(_apps):
    return settings.USER_AWARDS_ENABLED


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0231_user_awards"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["update_user_awards"]),
                        note="Backfill computed user titles and badges (#1819)",
                        test=_user_awards_enabled),
    ]
