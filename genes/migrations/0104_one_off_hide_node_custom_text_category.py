# The initial-data migration created NodeCustomText with hidden=False, so the hidden=True the node (#2386)
# code passes to get_or_create_category never applied and Node_*_custom lists leaked into gene list grids.
from django.db import migrations


def hide_node_custom_text(apps, _schema_editor):
    GeneListCategory = apps.get_model("genes", "GeneListCategory")
    GeneListCategory.objects.filter(name="NodeCustomText").update(hidden=True)


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0103_cached_web_resource_to_snpdb'),
    ]

    operations = [
        migrations.RunPython(hide_node_custom_text, migrations.RunPython.noop),
    ]
