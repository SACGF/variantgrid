from django.conf import settings
from django.db import migrations


def _set_global_default(apps, _schema):
    GlobalSettings = apps.get_model("snpdb", "GlobalSettings")
    global_settings = GlobalSettings.objects.get()
    if global_settings.node_grid_auto_load_max_variants is None:
        global_settings.node_grid_auto_load_max_variants = settings.ANALYSIS_NODE_GRID_AUTO_LOAD_MAX_VARIANTS
        global_settings.save()


def _reverse(apps, _schema):
    GlobalSettings = apps.get_model("snpdb", "GlobalSettings")
    global_settings = GlobalSettings.objects.get()
    global_settings.node_grid_auto_load_max_variants = None
    global_settings.save()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0190_settingsoverride_node_grid_auto_load_max_variants'),
    ]

    operations = [
        migrations.RunPython(_set_global_default, _reverse)
    ]
