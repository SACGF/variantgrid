from django.db import migrations
from django.db.models import Max

from annotation import vep_columns
from annotation.enums import VEPPlugin
from manual.operations.manual_operations import ManualOperation

# VariantAnnotationVersion pin field -> the plugin whose data file it pins
_PLUGIN_PINS = {
    "utr_annotator": VEPPlugin.UTR_ANNOTATOR,
    "avi": VEPPlugin.AVI,
    "splicevault": VEPPlugin.SPLICEVAULT,
}
BACKFILL_FLAG = "backfilled_utr_avi_splicevault"


def _versions_to_backfill(apps):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(**{BACKFILL_FLAG: False})


def _has_versions_to_backfill(apps) -> bool:
    return _versions_to_backfill(apps).exists()


def _backfill_plugin_columns(vav) -> dict[VEPPlugin, dict[str, str]]:
    """ VariantGrid column -> CSQ field, for each plugin this version pins whose columns it covers """
    plugin_columns = {}
    for field, vep_plugin in _PLUGIN_PINS.items():
        if getattr(vav, field):
            column_defs = vep_columns.filter_for(vep_plugin=vep_plugin, genome_build_name=vav.genome_build_id,
                                                 columns_version=vav.columns_version, vep_version=vav.vep)
            if column_defs:
                plugin_columns[vep_plugin] = {c.variant_grid_columns[0]: c.vep_info_field for c in column_defs}
    return plugin_columns


def _backfill_note(apps) -> str:
    """ The commands for each version, bounded by the last variant it had locked for annotation when the
        plugins were switched on - anything after that is annotated with them """
    AnnotationRangeLock = apps.get_model("annotation", "AnnotationRangeLock")
    lines = ["Fill in UTRAnnotator, AVI and SpliceVault annotation of variants annotated before they were switched "
             "on (#2426), with one VEP run per version. Full recipe, including downloading AVI.pm and the data "
             "files and converting SpliceVault with 'manage.py splicevault_to_refseq', in "
             "claude/research/annotation.md, 'A new column on an existing VariantAnnotationVersion'. "
             "For each version:"]
    for vav in _versions_to_backfill(apps).order_by("genome_build_id", "pk"):
        build = vav.genome_build_id
        plugin_columns = _backfill_plugin_columns(vav)
        columns = {column: field for c in plugin_columns.values() for column, field in c.items()}
        max_variant_id = AnnotationRangeLock.objects.filter(version=vav).aggregate(Max("max_variant_id"))
        common = f"--genome-build {build} --variant-annotation-version {vav.pk} --columns {','.join(columns)}"
        csq_args = " ".join(f"--csq {column}={field}" for column, field in columns.items())
        csq_fields = ",".join(["Feature", "PICK", *columns.values()])
        plugins = ", ".join(vep_plugin.label for vep_plugin in plugin_columns)
        lines.extend([
            "",
            f"{build} (VariantAnnotationVersion {vav.pk}):",
            f"1. python3 manage.py annotation_backfill_columns --dump {common} "
            f"--max-variant-id {max_variant_id['max_variant_id__max']} --output vep_plugins_{build}.vcf.gz",
            f"2. Annotate it with VEP: this build's STANDARD 'Show VEP command' on /annotation/annotation_versions, with "
            f"-i vep_plugins_{build}.vcf.gz -o vep_plugins_{build}.vep_annotated.vcf.gz, every --plugin and --custom "
            f"except {plugins} removed, and --fields {csq_fields}. Keep the other flags, so the PICK'd "
            f"transcript matches",
            f"3. python3 manage.py annotation_backfill_columns --import {common} {csq_args} "
            f"--vcf vep_plugins_{build}.vep_annotated.vcf.gz --dry-run",
            f"4. The same without --dry-run, plus --mark-backfilled {BACKFILL_FLAG}",
        ])
    return "\n".join(lines)


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0004_utr_annotator_avi_splicevault"),
        ("manual", "0001_initial"),
    ]

    operations = [
        ManualOperation.operation_other(args=["Backfill UTRAnnotator, AVI and SpliceVault columns by hand (#2426) "
                                              "- see note"],
                                        note=_backfill_note, test=_has_versions_to_backfill),
    ]
