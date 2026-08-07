from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0194_alter_vcflengthstats_variant_class'),
    ]

    operations = [
        migrations.DeleteModel(
            name='ColumnVCFInfo',
        ),
    ]
