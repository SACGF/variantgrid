# Generated by Django 3.1.3 on 2021-04-29 05:58

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('analysis', '0036_analysisnode_cloned_from'),
    ]

    operations = [
        migrations.AlterField(
            model_name='damagenode',
            name='impact_min',
            field=models.CharField(blank=True, choices=[('O', 'MODIFIER'), ('L', 'LOW'), ('M', 'MODERATE'), ('*', 'MODERATE*'), ('H', 'HIGH')], max_length=1, null=True),
        ),
    ]
