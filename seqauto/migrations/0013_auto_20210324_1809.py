# Generated by Django 3.1.3 on 2021-03-24 07:39

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('seqauto', '0012_auto_20210324_1702'),
    ]

    operations = [
        migrations.AddField(
            model_name='seqautomessage',
            name='code',
            field=models.TextField(null=True),
        ),
        migrations.AddField(
            model_name='seqautomessage',
            name='open',
            field=models.BooleanField(default=True),
        ),
        migrations.AlterField(
            model_name='seqautomessage',
            name='message',
            field=models.TextField(),
        ),
        migrations.AlterField(
            model_name='seqautomessage',
            name='record',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='seqauto.seqautorecord'),
        ),
        migrations.AlterField(
            model_name='seqautomessage',
            name='seqauto_run',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='seqauto.seqautorun'),
        ),
    ]
