from django.db import migrations


def add_flag_condition_resolution(apps, _schema_editor):
    flagtype = apps.get_model("flags", "FlagType")
    add_flag = flagtype(id='condition_resolution', label='Condition Resolution',
                        description='This classification has been recently modified.', context_id='classification',
                        comments_enabled=False, attributes={}, raise_permission='A', permission='A', importance=1)
    add_flag.save()


def add_flag_type_condition_resolution(apps, _schema_editor):
    flagtype_res = apps.get_model("flags", "FlagTypeResolution")
    add_flagtype_res = flagtype_res(flag_type_id='condition_resolution', resolution_id='closed')
    add_flagtype_res.save()


class Migration(migrations.Migration):

    dependencies = [
        ('flags', '0008_one_off_pending_change_flag'),
    ]

    operations = [
        migrations.RunPython(add_flag_condition_resolution),
        migrations.RunPython(add_flag_type_condition_resolution)
    ]
