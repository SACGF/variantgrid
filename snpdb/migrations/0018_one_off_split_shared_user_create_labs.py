from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_holding_organization(apps):
    holding_group_names = [org_group_name for org_group_name, lab_pattern in settings.USER_CREATE_ORG_LABS.items()
                           if "%(" in lab_pattern]
    return apps.get_model("snpdb", "Organization").objects.filter(group_name__in=holding_group_names).exists()


class Migration(migrations.Migration):

    dependencies = [
        ('manual', '0001_initial'),
        ('snpdb', '0017_lab_request'),
    ]

    operations = [
        ManualOperation.operation_manage(["one_off_split_shared_user_create_labs"],
                                         test=_test_has_holding_organization),
    ]
