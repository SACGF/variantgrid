# Generated by Django 4.1.3 on 2022-12-15 07:35

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_old_cdot(apps):
    TranscriptVersion = apps.get_model("genes", "TranscriptVersion")

    if not TranscriptVersion.objects.exists():
        return False  # No transcripts, mew install

    if cdot_transcript := TranscriptVersion.objects.filter(data__cdot__isnull=False).first():
        cdot_version = tuple(int(i) for i in cdot_transcript.data["cdot"].split("."))
        return cdot_version < (0, 2, 12)  # This is release with MANE/RefSeq etc tags
    return False


class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0061_remove_pfamsequenceidentifier_transcript_version'),
    ]

    operations = [
        ManualOperation.operation_other(args=["Update cdot - import_gene_annotation with cdot transcript DATA >= 0.2.12"
                                              " You can do this by running: "
                                              "./annotation/annotation_data/cdot_update.sh or manually, see"
                                              "https://github.com/SACGF/cdot/wiki/Download-JSON.gz-files"
                                              "Any cdot 0.2.X code is ok"],
                                        test=_test_old_cdot),
    ]