# Generated by Django 4.2.11 on 2024-08-07 03:48

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('classification', '0153_one_off_remind_rematch_allele_info_bad_end_norm'),
    ]

    operations = [
        migrations.AlterField(
            model_name='importedalleleinfo',
            name='status',
            field=models.CharField(choices=[('P', 'Processing'), ('I', 'Matched Imported Variant'), ('M', 'Complete'), ('F', 'Failed')], default='P', max_length=1),
        ),
    ]