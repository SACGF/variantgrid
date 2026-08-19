from django.db import migrations, models

# Issue #343 - tags in grid columns were unordered, making e.g. "artifact" hard to spot among many tags.
# Tags without an explicit sort_order sort as 0 (ties broken by tag name), so a large value puts artifact last
# by default while still letting users drag it elsewhere on the tag colors page.
ARTIFACT_SORT_ORDER = 1000


def _artifact_tag_sorts_last(apps, _schema_editor):
    Tag = apps.get_model("snpdb", "Tag")
    TagColor = apps.get_model("snpdb", "TagColor")
    TagColorsCollection = apps.get_model("snpdb", "TagColorsCollection")

    if not Tag.objects.filter(pk="artifact").exists():
        return

    for collection in TagColorsCollection.objects.all():
        tag_color, created = TagColor.objects.get_or_create(collection=collection, tag_id="artifact",
                                                            defaults={"rgb": "", "sort_order": ARTIFACT_SORT_ORDER})
        if not created and tag_color.sort_order is None:
            tag_color.sort_order = ARTIFACT_SORT_ORDER
            tag_color.save()
        collection.version_id += 1
        collection.save()


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0211_settingsoverride_variant_tag_stale_days'),
    ]

    operations = [
        migrations.AddField(
            model_name='tagcolor',
            name='sort_order',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.RunPython(_artifact_tag_sorts_last, migrations.RunPython.noop),
    ]
