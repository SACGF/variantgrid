"""
The cytoband files only hold the finest level (11p15.5 but no 11p15, 11p1 or 11p), so this generates a Cytoband
for every level above (#2685), spanning the rows underneath it with a blank stain, and records each row's level.
ISCN names are hierarchical: dropping the last digit of a name (and a trailing '.') gives its parent.
"""
from collections import defaultdict

from django.db import migrations, models

ARM = "A"
REGION = "R"
BAND = "B"
SUB_BAND = "S"
SUB_SUB_BAND = "SS"


def _cytoband_level(name: str) -> str:
    arm_and_band, _, sub_band = name.partition(".")
    if sub_band:
        return SUB_BAND if len(sub_band) == 1 else SUB_SUB_BAND
    return {1: ARM, 2: REGION}.get(len(arm_and_band), BAND)


def _parent_name(name: str) -> str:
    """ '' for an arm """
    return name[:-1].rstrip(".")


def _generate_cytoband_levels(apps, _schema_editor):
    Cytoband = apps.get_model("snpdb", "Cytoband")

    cytobands = list(Cytoband.objects.all())
    for cytoband in cytobands:
        cytoband.level = _cytoband_level(cytoband.name)
    Cytoband.objects.bulk_update(cytobands, ["level"], batch_size=2000)

    existing = {(c.contig_id, c.name) for c in cytobands}
    extents = defaultdict(list)  # (contig_id, name) -> [start, end]
    for cytoband in cytobands:
        name = cytoband.name
        while name := _parent_name(name):
            key = (cytoband.contig_id, name)
            if key in existing:
                continue
            if extent := extents.get(key):
                extent[0] = min(extent[0], cytoband.start)
                extent[1] = max(extent[1], cytoband.end)
            else:
                extents[key] = [cytoband.start, cytoband.end]

    Cytoband.objects.bulk_create([Cytoband(contig_id=contig_id, name=name, level=_cytoband_level(name),
                                           start=start, end=end, stain="")
                                  for (contig_id, name), (start, end) in extents.items()],
                                 batch_size=2000)


def _delete_generated_cytobands(apps, _schema_editor):
    Cytoband = apps.get_model("snpdb", "Cytoband")
    Cytoband.objects.filter(stain="").delete()


class Migration(migrations.Migration):

    dependencies = [
        ("snpdb", "0011_cytoband"),
    ]

    operations = [
        migrations.AddField(
            model_name="cytoband",
            name="level",
            field=models.CharField(
                choices=[
                    ("A", "Arm"),
                    ("R", "Region"),
                    ("B", "Band"),
                    ("S", "Sub-band"),
                    ("SS", "Sub-sub-band"),
                ],
                default="B",
                max_length=2,
            ),
            preserve_default=False,
        ),
        migrations.RunPython(_generate_cytoband_levels, reverse_code=_delete_generated_cytobands),
    ]
