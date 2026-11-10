from django.conf import settings
from django.db import migrations

from annotation.annotsv_regulatory_elements import (
    sort_instructions,
    unsorted_regulatory_element_files,
)
from manual.operations.manual_operations import ManualOperation


def _has_unsorted_regulatory_elements(_apps) -> bool:
    return settings.ANNOTATION_ANNOTSV_ENABLED and bool(unsorted_regulatory_element_files())


def _sort_note(_apps) -> str:
    return sort_instructions(unsorted_regulatory_element_files())


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0015_annotation_run_max_variant_length_tool_peak_memory"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        ManualOperation.operation_other(args=["Sort AnnotSV's regulatory element files by hand (#2749) - see note"],
                                        note=_sort_note, test=_has_unsorted_regulatory_elements),
    ]
