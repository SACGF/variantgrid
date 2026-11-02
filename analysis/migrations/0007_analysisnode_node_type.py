from django.db import migrations, models


def _backfill_node_type(apps, _schema_editor):
    """ Each node's concrete class name, from whichever subclass table holds its row. Deepest subclasses last, so a
        grandchild's name wins over its parent's. A row in no subclass table loads as a bare AnalysisNode """
    AnalysisNode = apps.get_model("analysis", "AnalysisNode")
    subclasses = [model for model in apps.get_app_config("analysis").get_models()
                  if AnalysisNode in model._meta.get_parent_list()]
    for model in sorted(subclasses, key=lambda m: len(m._meta.get_parent_list())):
        AnalysisNode.objects.filter(pk__in=model.objects.values("pk")).update(node_type=model.__name__)
    AnalysisNode.objects.filter(node_type__isnull=True).update(node_type="AnalysisNode")


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0006_damage_node_avi_splicevault_utr'),
    ]

    operations = [
        migrations.AddField(
            model_name='analysisnode',
            name='node_type',
            field=models.TextField(editable=False, null=True),
        ),
        migrations.RunPython(_backfill_node_type, migrations.RunPython.noop),
    ]
