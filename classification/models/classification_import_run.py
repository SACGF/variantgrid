from datetime import timedelta
from typing import Optional

from django.db.models.signals import post_save
from django import dispatch
from django.dispatch import receiver
from django.utils.timezone import now
from model_utils.models import TimeStampedModel
from django.db import models

from library.log_utils import NotificationBuilder

classification_imports_complete_signal = dispatch.Signal()


class ClassificationImportRunStatus(models.TextChoices):
    ONGOING = "O", 'Ongoing'
    COMPLETED = "C", 'Completed'
    UNFINISHED = "U", 'Unfinished'


class ClassificationImportRun(TimeStampedModel):
    # Could put lab and user in this to make it more specific
    # but an import can even be across labs
    identifier = models.TextField()
    row_count = models.IntegerField(default=0)
    status = models.TextField(choices=ClassificationImportRunStatus.choices, default=ClassificationImportRunStatus.ONGOING)

    @staticmethod
    def record_classification_import(identifier: str, add_row_count: int = 0, is_complete: bool = False):
        """
        :param identifier: An identifier for the import - when importing more up to date versions of the same file, try to re-use the same identifier
        :param add_row_count: How many rows just got added
        :param is_complete: Is the import complete
        This method may trigger classification_imports_complete_signal
        """
        cir = ClassificationImportRun._get_or_create_for_identifier(identifier)
        cir.row_count += add_row_count
        if is_complete:
            cir.status = ClassificationImportRunStatus.COMPLETED
        cir.save()

    @staticmethod
    def _get_or_create_for_identifier(identifier: str):
        ClassificationImportRun.cleanup()

        if existing := ClassificationImportRun.objects.filter(identifier=identifier, status=ClassificationImportRunStatus.ONGOING).first():
            return existing
        else:
            nb = NotificationBuilder("Import started")
            nb.add_markdown(f":golfer: Import started {identifier}")
            nb.send()
            return ClassificationImportRun(identifier=identifier)

    @staticmethod
    def cleanup():
        too_old = now() - timedelta(minutes=5)
        for unfinished in ClassificationImportRun.objects.filter(status=ClassificationImportRunStatus.ONGOING, modified__lte=too_old):
            unfinished.status = ClassificationImportRunStatus.UNFINISHED
            unfinished.save()

    @staticmethod
    def ongoing_import() -> Optional[str]:
        # should this check to see if there are any abandoned imports
        if ongoing := ClassificationImportRun.objects.filter(status=ClassificationImportRunStatus.ONGOING).first():
            return ongoing.identifier or "ongoing-import"
        return None


@receiver(post_save, sender=ClassificationImportRun)
def outstanding_import_check(sender, instance: ClassificationImportRun, **kwargs):
    if instance.status != ClassificationImportRunStatus.ONGOING:
        ongoing_imports = bool(ClassificationImportRun.ongoing_import())

        nb = NotificationBuilder("Import started")
        emoji = ":golf:" if instance.status == ClassificationImportRunStatus.COMPLETED else ":skunk:"
        nb.add_markdown(f"{emoji} Import {instance.get_status_display()} {instance.identifier}{' ongoing imports.' if ongoing_imports else ''}")
        nb.send()
        if not ongoing_imports:
            classification_imports_complete_signal.send(sender=ClassificationImportRun)
