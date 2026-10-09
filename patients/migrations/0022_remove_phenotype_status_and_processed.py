""" PhenotypeDescription.status was never set and DescriptionProcessingStatus had no columns. TextPhenotype.processed
    is "match_version is not null" since #2131: a sentence processed but never stamped (matched before #2131) is now
    awaiting a rematch, which the #2131 `--stale` ManualOperation (annotation 0191) asks for anyway (#2135) """
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ("patients", "0021_phenotype_description_owner"),
    ]

    operations = [
        migrations.RemoveField(model_name="phenotypedescription", name="status"),
        migrations.DeleteModel(name="DescriptionProcessingStatus"),
        migrations.RemoveField(model_name="textphenotype", name="processed"),
    ]
