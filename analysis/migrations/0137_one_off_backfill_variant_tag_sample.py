from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _has_variant_tags_needing_sample(apps):
    VariantTag = apps.get_model("analysis", "VariantTag")
    return VariantTag.objects.filter(sample__isnull=True, analysis__isnull=False).exists()


class Migration(migrations.Migration):
    dependencies = [
        ("analysis", "0136_varianttag_sample"),
        ("snpdb", "0251_one_off_tags_requiring_classification"),
    ]

    operations = [
        ManualOperation(task_id=ManualOperation.task_id_manage(["one_off_backfill_variant_tag_sample"]),
                        note="Work out which sample existing taggings were about, so they show up in the "
                             "Classify & Report tab (variantgrid_sapath#246)",
                        test=_has_variant_tags_needing_sample),
    ]
