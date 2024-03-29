# Generated by Django 3.1 on 2020-09-29 05:30

import django.utils.timezone
import model_utils.fields
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='EmailLog',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created', model_utils.fields.AutoCreatedField(default=django.utils.timezone.now, editable=False, verbose_name='created')),
                ('modified', model_utils.fields.AutoLastModifiedField(default=django.utils.timezone.now, editable=False, verbose_name='modified')),
                ('subject', models.TextField()),
                ('html', models.TextField()),
                ('text', models.TextField()),
                ('from_email', models.TextField(blank=True, null=True)),
                ('recipient_list', models.TextField()),
                ('probably_sent', models.BooleanField()),
                ('single_email', models.BooleanField()),
            ],
            options={
                'abstract': False,
            },
        ),
    ]
