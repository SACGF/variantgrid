from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_case_colliding_tags(apps):
    """ Tag names are the primary key, so 'artefact' and 'Artefact' are separate tags that split
        reporting and node filters - @see https://github.com/SACGF/variantgrid/issues/1751 """
    Tag = apps.get_model("snpdb", "Tag")
    lower_names = set()
    for tag_id in Tag.objects.values_list("pk", flat=True):
        lower_name = tag_id.lower()
        if lower_name in lower_names:
            return True
        lower_names.add(lower_name)
    return False


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0207_gene_level_contig"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["merge_tags", "--case-collisions",
                                                                "--delete-duplicate-variant-tags"]),
                        note="Merge tags that differ only by case (eg 'artefact' into 'Artefact') and delete "
                             "variant tags repeating the same variant/tag/analysis/user (#1751)",
                        test=_has_case_colliding_tags),
    ]
