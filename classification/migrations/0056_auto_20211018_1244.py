# Generated by Django 3.2.6 on 2021-10-18 02:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0055_new_ekey_lrg_identifier'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='clinvarexport',
            name='release_status',
        ),
        migrations.AlterField(
            model_name='clinvarexport',
            name='status',
            field=models.CharField(choices=[('N', 'New Submission'), ('C', 'Changes Pending'), ('D', 'Up to Date'), ('E', 'Error'), ('X', 'Exclude')], default='N', max_length=1),
        ),
    ]
