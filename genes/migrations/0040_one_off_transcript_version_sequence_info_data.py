# Generated by Django 3.2.6 on 2021-09-14 12:08

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_has_transcript_versions(apps):
    """ Don't need to run this for new deployments """

    TranscriptVersion = apps.get_model("genes", "TranscriptVersion")
    return TranscriptVersion.objects.exists()

class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0039_transcriptversionsequenceinfo_transcriptversionsequenceinfofastafileimport'),
    ]

    operations = [
        ManualOperation.operation_other(args=[
            "*** BEFORE rematching - import_transcript_fasta - see annotation page"],
            test=_test_has_transcript_versions),
    ]
