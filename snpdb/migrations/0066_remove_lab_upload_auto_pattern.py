# Generated by Django 4.0.2 on 2022-03-31 04:48

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0065_lab_contact_email_lab_contact_name_lab_contact_phone'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='lab',
            name='upload_auto_pattern',
        ),
    ]