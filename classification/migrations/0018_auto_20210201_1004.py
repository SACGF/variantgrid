# Generated by Django 3.1 on 2021-01-31 23:34

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0017_discordancereport_resolved_text'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='conditiontext',
            name='requires_approval',
        ),
        migrations.AddField(
            model_name='conditiontext',
            name='status',
            field=models.CharField(choices=[('T', 'Terms Provided'), ('A', 'Auto-Matched'), ('N', 'Not Auto-Matched'), ('U', 'User Reviewed')], default='N', max_length=1),
        ),
    ]