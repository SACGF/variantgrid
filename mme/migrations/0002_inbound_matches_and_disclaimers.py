import django.db.models.deletion
import django.utils.timezone
from django.db import migrations, models


class Migration(migrations.Migration):
    """ Per-peer attribution of inbound queries, the records we returned to each peer, and
        the disclaimer/terms a node returned with a response. No backfill - MME is off on
        every deployment, so these tables are empty. """

    dependencies = [
        ('classification', '0166_reminder_legacy_somatic_forwardport'),
        ('mme', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='mmeinboundquery',
            name='peer_node_id',
            field=models.CharField(blank=True, default='', max_length=64),
        ),
        migrations.AddField(
            model_name='mmesubmission',
            name='response_disclaimer',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.AddField(
            model_name='mmesubmission',
            name='response_terms',
            field=models.TextField(blank=True, default=''),
        ),
        migrations.CreateModel(
            name='MMEInboundMatch',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('score', models.FloatField()),
                ('remote_patient_id', models.CharField(blank=True, default='', max_length=255)),
                ('query_patient_json', models.JSONField()),
                ('notified', models.DateTimeField(blank=True, null=True)),
                ('created', models.DateTimeField(default=django.utils.timezone.now)),
                ('classification', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='classification.classification')),
                ('inbound_query', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='mme.mmeinboundquery')),
            ],
        ),
        migrations.AddIndex(
            model_name='mmeinboundmatch',
            index=models.Index(fields=['classification', 'remote_patient_id'],
                               name='mme_mmeinbo_classif_d46eaa_idx'),
        ),
    ]
