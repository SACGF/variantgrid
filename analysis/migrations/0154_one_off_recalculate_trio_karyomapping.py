from django.db import migrations


class Migration(migrations.Migration):
    """ Trio karyomapping is recounted lazily by analysis/0157 (#1754) """
    dependencies = [
        ("analysis", "0153_samplenode_copy_ratio_thresholds"),
    ]

    operations = []
