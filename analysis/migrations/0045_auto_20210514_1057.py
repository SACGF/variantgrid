# Generated by Django 3.2.1 on 2021-05-14 01:27

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0044_alter_zygositynode_zygosity'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='analysisnode',
            name='db_pid',
        ),
        migrations.AddField(
            model_name='nodetask',
            name='db_pid',
            field=models.IntegerField(null=True),
        ),
    ]
