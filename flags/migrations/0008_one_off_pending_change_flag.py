# Generated by Django 4.0.4 on 2022-06-23 05:43

from django.db import migrations


def _pending_changes_flag(apps, _schema_editor):
    FlagType = apps.get_model("flags", "FlagType")
    FlagResolution = apps.get_model("flags", "FlagResolution")
    FlagTypeResolution = apps.get_model("flags", "FlagTypeResolution")

    pending_changes_flag, _ = FlagType.objects.get_or_create(
        id="classification_pending_changes",
        defaults={
            "comments_enabled": False,
            "context_id": "classification",
            "description": "This classification has pending changes not yet reflected here",
            "importance": 2,
            "label": "Pending Changes",
            "only_one": True,
            "permission": "X",
            "raise_permission": "A"
        }
    )

    pending_changes, _ = FlagResolution.objects.get_or_create(
        id="pc_pending_changes",
        defaults={
            "label": "Pending Changes",
            "description": "This classification has pending changes not yet reflected here",
            "status": "O"
        }
    )
    changes_applied, _ = FlagResolution.objects.get_or_create(
        id="pc_applied",
        defaults={
            "label": "Changes Applied",
            "description": "Pending changes have been applied",
            "status": "C"
        }
    )
    changes_rejected, _ = FlagResolution.objects.get_or_create(
        id="pc_rejected",
        defaults={
            "label": "Changes Rejected",
            "description": "Pending changes were rejected",
            "status": "R"
        }
    )

    FlagTypeResolution.objects.get_or_create(
        flag_type=pending_changes_flag,
        resolution=pending_changes
    )

    FlagTypeResolution.objects.get_or_create(
        flag_type=pending_changes_flag,
        resolution=changes_applied
    )

    FlagTypeResolution.objects.get_or_create(
        flag_type=pending_changes_flag,
        resolution=changes_rejected
    )


class Migration(migrations.Migration):

    dependencies = [
        ('flags', '0007_alter_flagtype_permission_and_more'),
    ]

    operations = [
        migrations.RunPython(_pending_changes_flag)
    ]