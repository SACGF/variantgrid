from django.db import migrations
from django.db.models import Count

from manual.operations.manual_operations import ManualOperation

DUPLICATE_FIELDS = ["variant_id", "tag_id", "analysis_id", "user_id"]


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


def _has_duplicate_variant_tags(apps):
    """ Same user tagging the same variant with the same tag in the same analysis more than once """
    VariantTag = apps.get_model("analysis", "VariantTag")
    return VariantTag.objects.values(*DUPLICATE_FIELDS).annotate(count=Count("pk")).filter(count__gt=1).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0207_gene_level_contig"),
        ("analysis", "0112_tagnode_node_input"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["variant_tags", "merge-case-collisions"]),
                        note="Merge tags that differ only by case, eg 'artefact' into 'Artefact' (#1751)",
                        test=_has_case_colliding_tags),
        ManualOperation(task_id=ManualOperation.task_id_manage(["variant_tags", "delete-duplicates"]),
                        note="Delete variant tags repeating the same variant/tag/analysis/user (#1751)",
                        test=_has_duplicate_variant_tags),
    ]
