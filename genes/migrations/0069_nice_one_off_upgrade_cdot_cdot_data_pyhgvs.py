# Generated by Django 4.1.3 on 2023-05-18 07:14
from importlib import metadata
from importlib.metadata import PackageNotFoundError

from django.db import migrations

from manual.operations.manual_operations import ManualOperation


def _test_old_cdot_version(apps):
    needs_upgrade = False
    try:
        required_version = (0, 2, 17)  # Need FastaSeqFetcher fixes
        versions = [int(v) for v in metadata.version("cdot").split(".")]
        for req_v, v in zip(versions, required_version):
            needs_upgrade |= req_v < v
    except PackageNotFoundError:
        needs_upgrade = True

    return needs_upgrade


def _test_old_cdot_data(apps):
    TranscriptVersion = apps.get_model("genes", "TranscriptVersion")

    if not TranscriptVersion.objects.exists():
        return False  # No transcripts, new install

    if cdot_transcript := TranscriptVersion.objects.filter(data__cdot__isnull=False).first():
        cdot_version = tuple(int(i) for i in cdot_transcript.data["cdot"].split("."))
        return cdot_version < (0, 2, 17)  # 0.2.16 had chrom names instead of contigs
    return True  # ancient version of transcript annotations (no cdot version entry)


def _test_old_pyhgvs(apps):
    needs_upgrade = False
    try:
        required_version = (0, 12, 4)  # Bumped version on 12 Aug 2023 - allow spans to have = (match ref)
        versions = [int(v) for v in metadata.version("pyhgvs").split(".")]
        for req_v, v in zip(versions, required_version):
            needs_upgrade |= req_v < v
    except PackageNotFoundError:
        needs_upgrade = True

    return needs_upgrade



class Migration(migrations.Migration):

    dependencies = [
        ('genes', '0068_one_off_upgrade_cdot_data'),
    ]

    operations = [
        ManualOperation.operation_other(args=[
            "Update cdot library - sudo python3 -m pip install --upgrade cdot"],
            test=_test_old_cdot_version),
        ManualOperation.operation_other(args=["Update cdot - import_gene_annotation with cdot transcript DATA >= 0.2.17"
                                              " You can do this by running: "
                                              "./annotation/annotation_data/cdot_update.sh or manually, see "
                                              "https://github.com/SACGF/cdot/releases"
                                              "Data only - no need to update the library code for this change"],
                                        test=_test_old_cdot_data),
        ManualOperation.operation_other(args=[
            "Update pyHGVS library - sudo python3 -m pip install --force --upgrade git+https://github.com/SACGF/hgvs#egg=pyhgvs"],
            test=_test_old_pyhgvs),
    ]