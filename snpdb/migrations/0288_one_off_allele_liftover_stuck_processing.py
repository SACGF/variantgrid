from django.db import migrations


def _one_off_allele_liftover_stuck_processing(apps, _schema_editor):
    """ AlleleLiftovers left PROCESSING by runs that had finished - #2246
        - pre-liftover check failures were set back to PROCESSING when the run's VCF was made
        - a run whose upload pipeline failed or timed out after making its VCF never finished its alleles """
    AlleleLiftover = apps.get_model("snpdb", "AlleleLiftover")
    processing_qs = AlleleLiftover.objects.filter(status='P')
    processing_qs.filter(error__message__startswith="Validation Error").update(status='E')

    failed_pipeline_statuses = ['E', 'Z']  # ProcessingStatus ERROR, TIMED_OUT
    failed_run_qs = processing_qs.filter(
        liftover__uploadedliftover__file_upload__uploadpipeline__status__in=failed_pipeline_statuses)
    failed_run_qs.update(status='E', error={"message": "Liftover run failed"})


class Migration(migrations.Migration):

    dependencies = [
        ('snpdb', '0287_somalier_ancestry_batch'),
        ('upload', '0057_uploadstep_lease_expires'),
    ]

    operations = [
        migrations.RunPython(_one_off_allele_liftover_stuck_processing, migrations.RunPython.noop)
    ]
