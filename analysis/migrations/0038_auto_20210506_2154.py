# Generated by Django 3.2.1 on 2021-05-06 12:24

import django.db.models.deletion
import django_extensions.db.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0037_auto_20210429_1528'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='analysisnode',
            name='analysis_update_uuid',
        ),
        migrations.RemoveField(
            model_name='analysisnode',
            name='celery_task',
        ),
        migrations.CreateModel(
            name='NodeTask',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', django_extensions.db.fields.CreationDateTimeField(auto_now_add=True, verbose_name='created')),
                ('modified', django_extensions.db.fields.ModificationDateTimeField(auto_now=True, verbose_name='modified')),
                ('version', models.IntegerField()),
                ('analysis_update_uuid', models.UUIDField()),
                ('celery_task', models.CharField(max_length=36, null=True)),
                ('node', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='analysis.analysisnode')),
            ],
            options={
                'unique_together': {('node', 'version')},
            },
        ),
    ]
