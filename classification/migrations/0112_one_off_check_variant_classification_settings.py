# Generated by Django 4.2.2 on 2023-09-26 05:13
import logging

from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _check_settings_for_variant_classification(apps):
    has_old_settings = False
    for x in dir(settings):
        if x.startswith("VARIANT_CLASSIFICATION"):
            if not has_old_settings:  # Only show 1st time
                logging.warning("Your developer settings has obsolete settings. Please rename:")
                has_old_settings = True
            old_setting = x
            new_setting = x[len("VARIANT_"):]
            logging.warning("settings.%s -> %s", old_setting, new_setting)
    return has_old_settings


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0111_alter_resolvedvariantinfo_options'),
    ]

    operations = [
        ManualOperation.operation_other([
            "Rename settings from VARIANT_CLASSIFICATION_X to CLASSIFICATION_X"
        ], test=_check_settings_for_variant_classification),
    ]