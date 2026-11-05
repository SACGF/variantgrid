""" The boolean "require zygosity" fields become the three-way ZygosityRequirement choice (#2659). Saved
node versions keep their meaning: True was "call required", False was "no call allowed". """

from django.db import migrations, models

# ZygosityRequirement codes, as literals so the migration stays independent of the enum
CALL_REQUIRED = "C"
ALLOW_NO_CALL = "N"

# (model, old boolean field, new choice field)
FIELDS = [
    ("DuoNode", "require_zygosity", "relative_zygosity_requirement"),
    ("PedigreeNode", "require_zygosity", "zygosity_requirement"),
    ("QuadNode", "require_parent_zygosity", "parent_zygosity_requirement"),
    ("QuadNode", "require_sibling_zygosity", "sibling_zygosity_requirement"),
    ("TrioNode", "require_zygosity", "parent_zygosity_requirement"),
]


def booleans_to_requirement(apps, _schema_editor):
    for model_name, old_field, new_field in FIELDS:
        model = apps.get_model("analysis", model_name)
        model.objects.filter(**{old_field: False}).update(**{new_field: ALLOW_NO_CALL})


def requirement_to_booleans(apps, _schema_editor):
    for model_name, old_field, new_field in FIELDS:
        model = apps.get_model("analysis", model_name)
        model.objects.filter(**{new_field: ALLOW_NO_CALL}).update(**{old_field: False})


def _requirement_field():
    return models.CharField(
        choices=[("N", "No call allowed"), ("C", "Call required"), ("R", "Called REF required")],
        default=CALL_REQUIRED, max_length=1,
    )


class Migration(migrations.Migration):

    dependencies = [
        ("analysis", "0011_population_node_presets"),
    ]

    operations = [
        migrations.AddField(model_name=model_name.lower(), name=new_field, field=_requirement_field())
        for model_name, _old_field, new_field in FIELDS
    ] + [
        migrations.RunPython(booleans_to_requirement, requirement_to_booleans),
    ] + [
        migrations.RemoveField(model_name=model_name.lower(), name=old_field)
        for model_name, old_field, _new_field in FIELDS
    ]
