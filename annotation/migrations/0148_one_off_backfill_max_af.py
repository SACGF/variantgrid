from django.db import migrations


# No-op: this migration originally registered a ManualOperation to run the
# `fix_historical_max_af` management command (backfilling VariantAnnotation.max_af
# for the PopulationNode optimisation, #1547). That whole optimisation was
# reverted - the max_af field is removed again in 0152 and the management command
# no longer exists - so registering the manual task would point operators at a
# command that isn't there. The migration node is kept (as a dependency anchor)
# but does nothing.
class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0147_variantannotation_max_af"),
    ]

    operations = [
    ]
