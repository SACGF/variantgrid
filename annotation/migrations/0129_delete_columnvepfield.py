from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('annotation', '0128_genepubmedcount_delete_genesymbolpubmedcount'),
    ]

    operations = [
        migrations.DeleteModel(
            name='ColumnVEPField',
        ),
    ]
