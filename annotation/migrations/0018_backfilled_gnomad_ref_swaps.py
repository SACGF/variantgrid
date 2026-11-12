from django.conf import settings
from django.db import migrations, models

from annotation.vep_config import parse_gnomad_version_from_filename
from manual.operations.manual_operations import ManualOperation

T2T = "T2T-CHM13v2.0"
GNOMAD_LIFTOVER = "4.1.t2t_liftover"
GNOMAD_REF_SWAPS = "4.1.t2t_liftover_ref_swaps"


def _settings_gnomad_version():
    gnomad_filename = settings.ANNOTATION.get(T2T, {}).get("vep_config", {}).get("gnomad4")
    return gnomad_filename and parse_gnomad_version_from_filename(gnomad_filename)


def _move_active_t2t_version_to_ref_swaps(apps, _schema_editor):
    """ Settings now pass VEP the ref swaps file, so the version's gnomAD pin follows them (as VEP checks it) and
        the variants annotated from the old file wait for the backfill. The old file's AFs are all in the new one
        unchanged, so a common partition built from it is valid for both (#1582) """
    if _settings_gnomad_version() != GNOMAD_REF_SWAPS:
        return
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    VariantAnnotationVersion.objects.filter(genome_build_id=T2T, status="ACTIVE", gnomad=GNOMAD_LIFTOVER) \
        .update(gnomad=GNOMAD_REF_SWAPS, backfilled_gnomad_ref_swaps=False)

    CohortGenotypeCommonFilterVersion = apps.get_model("snpdb", "CohortGenotypeCommonFilterVersion")
    for common_filter in CohortGenotypeCommonFilterVersion.objects.filter(genome_build_id=T2T,
                                                                          gnomad_version=GNOMAD_LIFTOVER):
        if GNOMAD_REF_SWAPS not in common_filter.additional_gnomad_versions:
            common_filter.additional_gnomad_versions.append(GNOMAD_REF_SWAPS)
            common_filter.save(update_fields=["additional_gnomad_versions"])


def _has_version_to_backfill(apps) -> bool:
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(backfilled_gnomad_ref_swaps=False).exists()


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0017_reannotate_converted_symbolic_variants"),
        ("snpdb", "0015_settingsoverride_analysis_node_zoom"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        migrations.AddField(
            model_name="variantannotationversion",
            name="backfilled_gnomad_ref_swaps",
            field=models.BooleanField(default=True),
        ),
        migrations.RunPython(_move_active_t2t_version_to_ref_swaps, migrations.RunPython.noop),
        ManualOperation.operation_manage(["one_off_backfill_t2t_gnomad_ref_swaps"],
                                         note="Fills gnomAD on T2T variants annotated before the 4.1.t2t_liftover_ref_swaps "
                                              "file (#2788). Needs the new file (download_annotation_data) and "
                                              "bcftools; streams the 50 GB file once, about an hour",
                                         test=_has_version_to_backfill),
    ]
