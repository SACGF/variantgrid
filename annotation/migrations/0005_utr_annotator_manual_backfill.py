from django.db import migrations
from django.db.models import Max

from manual.operations.manual_operations import ManualOperation

_COLUMNS_TO_CSQ_FIELDS = {
    "utr_annotator_consequence": "5UTR_consequence",
    "utr_annotator_annotation": "5UTR_annotation",
    "utr_annotator_existing_uorfs": "Existing_uORFs",
    "utr_annotator_existing_out_of_frame_oorfs": "Existing_OutOfFrame_oORFs",
    "utr_annotator_existing_inframe_oorfs": "Existing_InFrame_oORFs",
}


def _versions_to_backfill(apps):
    VariantAnnotationVersion = apps.get_model("annotation", "VariantAnnotationVersion")
    return VariantAnnotationVersion.objects.filter(backfilled_utr_annotator=False)


def _has_versions_to_backfill(apps) -> bool:
    return _versions_to_backfill(apps).exists()


def _backfill_note(apps) -> str:
    """ The commands for each version, bounded by the last variant it had locked for annotation when the
        plugin was switched on - anything after that is annotated with it """
    AnnotationRangeLock = apps.get_model("annotation", "AnnotationRangeLock")
    columns = ",".join(_COLUMNS_TO_CSQ_FIELDS)
    csq_args = " ".join(f"--csq {column}={field}" for column, field in _COLUMNS_TO_CSQ_FIELDS.items())
    csq_fields = ",".join(["Feature", "PICK", *_COLUMNS_TO_CSQ_FIELDS.values()])

    lines = ["Fill in UTRAnnotator (5' UTR / uORF) annotation of variants annotated before it was switched on "
             "(#2426). Full recipe in claude/research/annotation.md, 'A new column on an existing "
             "VariantAnnotationVersion'. First download uORF_5UTR_<build>_PUBLIC.txt from "
             "https://github.com/Ensembl/UTRannotator to the vep_config 'utr_annotator' path. For each version:"]
    for vav in _versions_to_backfill(apps).order_by("genome_build_id", "pk"):
        build = vav.genome_build_id
        max_variant_id = AnnotationRangeLock.objects.filter(version=vav).aggregate(Max("max_variant_id"))
        common = f"--genome-build {build} --variant-annotation-version {vav.pk} --columns {columns}"
        lines.extend([
            "",
            f"{build} (VariantAnnotationVersion {vav.pk}):",
            f"1. python3 manage.py annotation_backfill_columns --dump {common} --consequence 5_prime_UTR_variant "
            f"--max-variant-id {max_variant_id['max_variant_id__max']} --output utr_annotator_{build}.vcf.gz",
            "2. Annotate it with VEP: the pipeline_command of one of this version's STANDARD AnnotationRuns, "
            f"with -i utr_annotator_{build}.vcf.gz -o utr_annotator_{build}.vep_annotated.vcf.gz, every --plugin "
            f"and --custom replaced by the UTRAnnotator one, and --fields {csq_fields}. Keep the other flags, "
            "so the PICK'd transcript matches",
            f"3. python3 manage.py annotation_backfill_columns --import {common} {csq_args} "
            f"--vcf utr_annotator_{build}.vep_annotated.vcf.gz --dry-run",
            "4. The same without --dry-run, plus --mark-backfilled backfilled_utr_annotator",
        ])
    return "\n".join(lines)


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0004_utr_annotator"),
        ("manual", "0001_initial"),
    ]

    operations = [
        ManualOperation.operation_other(args=["Backfill UTRAnnotator columns by hand (#2426) - see note"],
                                        note=_backfill_note, test=_has_versions_to_backfill),
    ]
