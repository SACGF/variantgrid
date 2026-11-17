from django.conf import settings
from django.db import migrations

from annotation.annotsv_annotation import (
    get_annotsv_fork_fixes_installed,
    get_annotsv_fork_install_instructions,
    is_annotsv_fork_installed,
)
from manual.operations.manual_operations import ManualOperation


def _runs_older_annotsv_fork(_apps) -> bool:
    """ Stock AnnotSV is 0020's step """
    return (settings.ANNOTATION_ANNOTSV_ENABLED and bool(get_annotsv_fork_fixes_installed())
            and not is_annotsv_fork_installed())


def _fork_note(_apps) -> str:
    queue = settings.ANNOTATION_ANNOTSV_WORKER_QUEUE
    return ("The AnnotSV fork's vg4 tag also stops the output step holding every split line in memory, which still "
            "killed ranges with very many large SVs, and fixes ACMG criteria 4O/2C/2F never being applied, which "
            "left deletions and duplications with the genes of a common benign CNV as VUS (#2841). That is a "
            "correctness fix: the ACMG classes already stored are wrong for those SVs, so the fork reports its own "
            "-version and is a new pipeline version, and promoting it re-annotates every range. Install, register "
            "and promote in one sitting: AnnotSV runs are refused in between.\n"
            f"{get_annotsv_fork_install_instructions()}\n"
            f"AnnotSV runs now go to the {queue} celery queue (ANNOTATION_ANNOTSV_WORKER_QUEUE): enable and start "
            f"celeryd_{queue} on the machine that has AnnotSV installed (systemctl enable --now celeryd_{queue}), "
            "or no AnnotSV run is dispatched.")


class Migration(migrations.Migration):

    dependencies = [
        ("annotation", "0020_annotsv_fork"),
    ]

    operations = [
        ManualOperation.operation_other(args=["Install the AnnotSV fork's vg4 tag, re-run AnnotSV and start its celery worker (#2841) - see note"],
                                        note=_fork_note, test=_runs_older_annotsv_fork),
    ]
