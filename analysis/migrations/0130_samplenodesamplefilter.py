from django.db import migrations, models

THRESHOLD_FIELDS = ["min_ad", "min_dp", "min_gq", "max_pl"]


def _inherit_node_thresholds(apps, _schema_editor):
    """ Rows were written with the node's value in every field it wasn't overriding - null those so
        they follow the node the way a row saved from here on does """
    SampleNodeSampleFilter = apps.get_model("analysis", "SampleNodeSampleFilter")
    for row in SampleNodeSampleFilter.objects.select_related("node"):
        modified = False
        for field in THRESHOLD_FIELDS:
            if getattr(row, field) == getattr(row.node, field):
                setattr(row, field, None)
                modified = True
        if modified:
            row.save()


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0129_analysistemplateversion_active_default_false'),
    ]

    operations = [
        migrations.RenameModel(
            old_name='SampleNodeSampleThreshold',
            new_name='SampleNodeSampleFilter',
        ),
        migrations.AlterField(
            model_name='samplenodesamplefilter',
            name='min_ad',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='samplenodesamplefilter',
            name='min_dp',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='samplenodesamplefilter',
            name='min_gq',
            field=models.IntegerField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='samplenodesamplefilter',
            name='zygosity',
            field=models.CharField(blank=True, max_length=4, null=True),
        ),
        migrations.AddField(
            model_name='samplenodesamplefilter',
            name='pass_only',
            field=models.BooleanField(null=True),
        ),
        migrations.AddField(
            model_name='samplenodesamplefilter',
            name='af_min',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AddField(
            model_name='samplenodesamplefilter',
            name='af_max',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.RunPython(_inherit_node_thresholds, migrations.RunPython.noop),
    ]
