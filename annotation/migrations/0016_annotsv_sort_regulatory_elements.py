from django.conf import settings
from django.db import migrations

from annotation.annotsv_annotation import get_annotsv_sort_regulatory_elements_script
from manual.operations.manual_operations import ManualOperation


def _has_annotsv_bundle_with_unsorted_files(_apps) -> bool:
    return settings.ANNOTATION_ANNOTSV_ENABLED and settings.ANNOTATION_ANNOTSV_BUNDLE_VERSION.startswith("3.5")


def _sort_note(_apps) -> str:
    return ("The AnnotSV 3.5 annotations bundle ships unsorted regulatory element files, which cost every AnnotSV "
            "run about 1.4 GB more memory (#2749). AnnotSV's output is the same once they are sorted. Run this as "
            "the owner of the annotation files, on the machine that holds them:\n"
            f"{get_annotsv_sort_regulatory_elements_script()}")


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0015_annotation_run_max_variant_length_tool_peak_memory"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        ManualOperation.operation_other(args=["Sort AnnotSV's regulatory element files by hand (#2749) - see note"],
                                        note=_sort_note, test=_has_annotsv_bundle_with_unsorted_files),
    ]
