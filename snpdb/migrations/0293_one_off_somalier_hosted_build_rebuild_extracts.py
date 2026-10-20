from django.conf import settings
from django.db import migrations

from manual.operations.manual_operations import ManualOperation

_DOWNLOAD_ARGS = ["download_annotation_data"]
_REBUILD_ARGS = ["somalier_existing_vcfs", "--clear"]


def _has_somalier_extracts(apps):
    """ Extracts written for v0.3.5 have the AD pair in the site's order, which the somalier we now
        host doesn't want - ancestry comes back AMR for everyone (#183) """
    SomalierVCFExtract = apps.get_model("snpdb", "SomalierVCFExtract")
    return settings.SOMALIER["enabled"] and SomalierVCFExtract.objects.exists()


class Migration(migrations.Migration):
    dependencies = [
        ("snpdb", "0292_phastcons_mammal_labels"),
    ]

    operations = [
        ManualOperation.operation_manage(args=_DOWNLOAD_ARGS,
                                         note="Fetch the somalier build in settings.SOMALIER['annotation']['command'] "
                                              "from variantgrid.com (#183)",
                                         test=_has_somalier_extracts),
        ManualOperation.operation_manage(args=_REBUILD_ARGS,
                                         note="Rebuild somalier extracts without the allele order compensation, "
                                              "then relate and ancestry (#183)",
                                         test=_has_somalier_extracts,
                                         requires=[f"after:{ManualOperation.task_id_manage(_DOWNLOAD_ARGS)}"]),
    ]
