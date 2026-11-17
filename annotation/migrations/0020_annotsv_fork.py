from django.conf import settings
from django.db import migrations

from annotation.annotsv_annotation import (
    get_annotsv_fork_fixes_installed,
    get_annotsv_fork_install_instructions,
)
from manual.operations.manual_operations import ManualOperation


def _runs_stock_annotsv(_apps) -> bool:
    return settings.ANNOTATION_ANNOTSV_ENABLED and not get_annotsv_fork_fixes_installed()


def _fork_note(_apps) -> str:
    return ("Stock AnnotSV's regulatory elements step is killed for memory on ranges with many large SVs, on every "
            "retry (#2841). Our fork fixes it.\n"
            f"{get_annotsv_fork_install_instructions()}\n"
            "Promoting the new pipeline version creates fresh runs for every range, including those that failed for "
            "memory. Retrying or deleting a failed run removes what its attempts left in "
            f"{settings.ANNOTATION_VCF_DUMP_DIR}/annotsv_<run>, often several GB.")


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0019_new_t2t_version_gnomad_ref_swaps"),
        ("manual", "0002_seed_data"),
    ]

    operations = [
        ManualOperation.operation_other(args=["Install AnnotSV from our fork (#2841) - see note"],
                                        note=_fork_note, test=_runs_stock_annotsv),
    ]
