from django.conf import settings
from django.db import migrations

from annotation.vep_config import parse_gnomad_version_from_filename
from manual.operations.manual_operations import ManualOperation

T2T = "T2T-CHM13v2.0"
GNOMAD_LIFTOVER = "4.1.t2t_liftover"
GNOMAD_REF_SWAPS = "4.1.t2t_liftover_ref_swaps"


def _settings_gnomad_version():
    gnomad_filename = settings.ANNOTATION.get(T2T, {}).get("vep_config", {}).get("gnomad4")
    return gnomad_filename and parse_gnomad_version_from_filename(gnomad_filename)


def _move_new_t2t_versions_to_ref_swaps(apps, _schema_editor):
    """ 0018 moved only the ACTIVE T2T version: a NEW one (annotated ahead of switching it on) kept the old gnomAD
        pin and fails VEP's version check now settings pass the ref swaps file (#2837) """
    if _settings_gnomad_version() != GNOMAD_REF_SWAPS:
        return
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    VariantAnnotationVersion.objects.filter(genome_build_id=T2T, status__in=["ACTIVE", "NEW"],
                                            gnomad=GNOMAD_LIFTOVER) \
        .update(gnomad=GNOMAD_REF_SWAPS, backfilled_gnomad_ref_swaps=False)


def _has_version_to_backfill(apps) -> bool:
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(backfilled_gnomad_ref_swaps=False).exists()


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0018_backfilled_gnomad_ref_swaps"),
    ]

    operations = [
        migrations.RunPython(_move_new_t2t_versions_to_ref_swaps, migrations.RunPython.noop),
        # Registered again: a deployment that already ran it for the ACTIVE version would otherwise skip the NEW one
        ManualOperation.operation_manage(["one_off_backfill_t2t_gnomad_ref_swaps"],
                                         note="Fills gnomAD on T2T variants annotated before the 4.1.t2t_liftover_ref_swaps "
                                              "file (#2788, #2837). Needs the new file (download_annotation_data) and "
                                              "bcftools; streams the 50 GB file once, about an hour",
                                         test=_has_version_to_backfill),
    ]
