from django.db import migrations


def _one_off_allele_liftover_precheck_error_status(apps, _schema_editor):
    """ Pre-liftover check failures were set back to PROCESSING when the run's VCF was made - #2246 """
    AlleleLiftover = apps.get_model("snpdb", "AlleleLiftover")
    AlleleLiftover.objects.filter(status='P',
                                  error__message__startswith="Validation Error").update(status='E')


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0287_somalier_ancestry_batch'),
    ]

    operations = [
        migrations.RunPython(_one_off_allele_liftover_precheck_error_status, migrations.RunPython.noop)
    ]
