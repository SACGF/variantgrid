import django.db.models.deletion
from django.db import migrations, models


def _clear_node_tasks(apps, schema_editor):
    """ NodeTask is a transient lease/lock row - safe to drop before re-keying onto NodeVersion.
        Any analysis in flight at deploy is re-dispatched by reschedule_stalled_analyses / a user
        edit (its nodes are still DIRTY/loading). """
    apps.get_model("analysis", "NodeTask").objects.all().delete()


class Migration(migrations.Migration):

    dependencies = [
        ("analysis", "0104_alter_quadnode_inheritance_and_more"),
    ]

    operations = [
        # Re-key NodeTask from (node FK + version) onto a OneToOne against NodeVersion so it
        # shares that row's lifecycle (cascade cleanup of stale versions) - issue #346.
        migrations.RunPython(_clear_node_tasks, migrations.RunPython.noop),
        migrations.AlterUniqueTogether(
            name="nodetask",
            unique_together=set(),
        ),
        # Dropping the local Meta restores TimeStampedModel's inherited option
        migrations.AlterModelOptions(
            name="nodetask",
            options={"get_latest_by": "modified"},
        ),
        migrations.RemoveField(
            model_name="nodetask",
            name="node",
        ),
        migrations.RemoveField(
            model_name="nodetask",
            name="version",
        ),
        migrations.AddField(
            model_name="nodetask",
            name="node_version",
            field=models.OneToOneField(
                on_delete=django.db.models.deletion.CASCADE,
                to="analysis.nodeversion",
            ),
        ),
        # Lease + backoff fields (mocha lease record)
        migrations.AddField(
            model_name="nodetask",
            name="leased_by",
            field=models.CharField(max_length=64, null=True),
        ),
        migrations.AddField(
            model_name="nodetask",
            name="lease_expires",
            field=models.DateTimeField(null=True),
        ),
        migrations.AddField(
            model_name="nodetask",
            name="attempt_count",
            field=models.IntegerField(default=0),
        ),
        migrations.AddField(
            model_name="nodetask",
            name="last_attempt",
            field=models.DateTimeField(null=True),
        ),
        migrations.AddField(
            model_name="nodetask",
            name="run_after",
            field=models.DateTimeField(null=True),
        ),
    ]
