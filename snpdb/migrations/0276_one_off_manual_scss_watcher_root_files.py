from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _is_dev_box(_apps):
    return settings.DEBUG


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0275_remove_combined_variant_output_vcf_source_settings"),
    ]

    operations = [
        ManualOperation.operation_other(args=[
            "PyCharm: Settings > Tools > File Watchers > SCSS - tick 'Track only root files' and make sure the watcher runs Dart Sass (libsass/sassc can't read @use). global.scss is now @use lines over partials in css/global/ - see https://github.com/SACGF/variantgrid/issues/1906",
        ], test=_is_dev_box),
    ]
