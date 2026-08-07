from django.db import migrations

# analysis_create_default_templates imports analyses, and creating an Analysis needs a valid
# AnnotationVersion - so on a deployment that hasn't imported annotation yet the task fails.
# The gate is declared at the call sites (analysis migrations 0068/0080), which only stamps
# deployments that apply them from here on - this covers task rows created before that.
TASK_ID = "manage*analysis_create_default_templates"
REQUIRES = ["annotation-version"]


def add_requires(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationTask.objects.filter(pk=TASK_ID).update(requires=REQUIRES)


def clear_requires(apps, schema_editor):
    ManualMigrationTask = apps.get_model("manual", "ManualMigrationTask")
    ManualMigrationTask.objects.filter(pk=TASK_ID).update(requires=list())


class Migration(migrations.Migration):

    dependencies = [
        # 0108 follows 0080, which creates the task - so the update can't run before the row exists
        ('analysis', '0108_intersectionnodecontig'),
        ('manual', '0005_backfill_manual_task_requires'),  # ManualMigrationTask.requires
    ]

    operations = [
        migrations.RunPython(add_requires, clear_requires),
    ]
