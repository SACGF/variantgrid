# Generated by Django 4.2.5 on 2023-11-30 03:17

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0106_useraward'),
    ]

    operations = [
        migrations.AddField(
            model_name='useraward',
            name='award_level',
            field=models.TextField(choices=[('G', 'Gold'), ('S', 'Silver'), ('B', 'Bronze')], default='G', max_length=1),
        ),
    ]