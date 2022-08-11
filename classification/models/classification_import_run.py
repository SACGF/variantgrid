from datetime import timedelta
from typing import Optional

from django import dispatch
from django.db import models
from django.db.models import CASCADE
from django.db.models.signals import post_save
from django.dispatch import receiver
from django.utils.timezone import now
from model_utils.models import TimeStampedModel
from classification.models.classification_utils import ClassificationPatchStatus
from classification.models.uploaded_classifications_unmapped import UploadedClassificationsUnmapped
from library.log_utils import NotificationBuilder

classification_imports_complete_signal = dispatch.Signal()


class ClassificationImportRunStatus(models.TextChoices):
    ONGOING = "O", 'Ongoing'
    COMPLETED = "C", 'Completed'
    UNFINISHED = "U", 'Unfinished'


MAX_IMPORT_AGE = timedelta(minutes=5)


class ClassificationImportRun(TimeStampedModel):
    # Could put lab and user in this to make it more specific
    # but an import can even be across labs
    identifier = models.TextField()
    row_count = models.IntegerField(default=0)
    status = models.TextField(choices=ClassificationImportRunStatus.choices, default=ClassificationImportRunStatus.ONGOING)
    from_file = models.ForeignKey(UploadedClassificationsUnmapped, on_delete=CASCADE, null=True, blank=True)

    row_count_new = models.IntegerField(default=0)
    row_count_update = models.IntegerField(default=0)
    row_count_no_change = models.IntegerField(default=0)
    row_count_withdrawn = models.IntegerField(default=0)
    row_count_delete = models.IntegerField(default=0)
    row_count_un_withdrawn = models.IntegerField(default=0)
    row_count_already_withdrawn = models.IntegerField(default=0)
    row_count_unknown = models.IntegerField(default=0)
    missing_row_count = models.IntegerField(null=True, blank=True)
    logging_version = models.IntegerField(default=0)

    def __str__(self):
        parts = []
        parts.append(f"(run_id={self.pk})")
        if file := self.from_file:
            parts.append(file.filename)
        else:
            parts.append(self.identifier)
        parts.append(f"(rows:{self.row_count})")
        return "".join(parts)

    def apply_missing_row_count(self) -> int:
        if from_file := self.from_file:
            from classification.models import Classification
            self.missing_row_count = Classification.objects.filter(lab=from_file.lab, withdrawn=False).exclude(last_import_run__from_file=from_file).count()

    def increment_status(self, status: ClassificationPatchStatus):
        self.row_count += 1
        if status == ClassificationPatchStatus.NEW:
            self.row_count_new += 1
        elif status == ClassificationPatchStatus.ALREADY_WITHDRAWN:
            self.row_count_already_withdrawn += 1
        elif status == ClassificationPatchStatus.UPDATE:
            self.row_count_update += 1
        elif status == ClassificationPatchStatus.DELETED:
            self.row_count_delete += 1
        elif status == ClassificationPatchStatus.NO_CHANGE:
            self.row_count_no_change += 1
        elif status == ClassificationPatchStatus.WITHDRAWN:
            self.row_count_withdrawn += 1
        else:
            # shouldn't happen but it's here if we need it
            self.row_count_unknown += 1

    @staticmethod
    def record_classification_import(identifier: str, add_row_count: int = 0, is_complete: bool = False) -> 'ClassificationImportRun':
        """
        :param identifier: An identifier for the import - when importing more up to date versions of the same file, try to re-use the same identifier
        :param add_row_count: How many rows just got added
        :param is_complete: Is the import complete
        This method may trigger classification_imports_complete_signal
        """
        use_me: Optional[ClassificationImportRun] = None

        # see if there's already an ongoing import, if it's not too old
        too_old = now() - MAX_IMPORT_AGE
        use_me = ClassificationImportRun.objects.filter(status=ClassificationImportRunStatus.ONGOING, identifier=identifier, modified__gte=too_old).order_by('-modified').first()

        if not use_me:
            use_me = ClassificationImportRun(identifier=identifier, logging_version=1)
            nb = NotificationBuilder("Import started")
            nb.add_markdown(f":golfer: Import Started {identifier}")
            nb.send()

        use_me.row_count += add_row_count
        if is_complete:
            use_me.status = ClassificationImportRunStatus.COMPLETED
        use_me.save()

        # after we've either updated or created a new import, cleanup any old ones
        # do this after so we don't close and import and start a new one - but trigger the fact that there were 0 imports in that split moment inbetween
        ClassificationImportRun.cleanup()
        return use_me

    @staticmethod
    def cleanup():
        too_old = now() - MAX_IMPORT_AGE
        for unfinished in ClassificationImportRun.objects.filter(status=ClassificationImportRunStatus.ONGOING, modified__lt=too_old):
            unfinished.status = ClassificationImportRunStatus.UNFINISHED
            unfinished.save()

    @staticmethod
    def ongoing_imports() -> Optional[str]:
        # should this check to see if there are any abandoned imports
        if ongoings := list(ClassificationImportRun.objects.filter(status=ClassificationImportRunStatus.ONGOING).order_by('-created')):
            return ", ".join(ongoing.identifier for ongoing in ongoings) or "ongoing-import"
        return None


@receiver(post_save, sender=ClassificationImportRun)
def outstanding_import_check(sender, instance: ClassificationImportRun, **kwargs):
    if instance.status != ClassificationImportRunStatus.ONGOING:
        ongoing_imports = ClassificationImportRun.ongoing_imports()

        nb = NotificationBuilder("Import started")
        emoji = ":golf:" if instance.status == ClassificationImportRunStatus.COMPLETED else ":skunk:"
        ongoing_message = f" ongoing imports {ongoing_imports}" if ongoing_imports else ""
        nb.add_markdown(f"{emoji} Import {instance.get_status_display()} {instance.identifier} {instance.row_count} rows{ongoing_message if ongoing_imports else ''}")
        # provide full details of import numbers in notification
        if instance.row_count_new:
            nb.add_field("New", instance.row_count_new)
        if instance.row_count_update:
            nb.add_field("Updated", instance.row_count_update)
        if instance.row_count_no_change:
            nb.add_field("No Change", instance.row_count_no_change)
        # less likely, but still report them just in case
        if instance.row_count_un_withdrawn:
            nb.add_field("Un-Withdrawn", instance.row_count_un_withdrawn)
        if instance.row_count_already_withdrawn:
            nb.add_field("Already Withdrawn", instance.row_count_already_withdrawn)
        if instance.row_count_withdrawn:
            nb.add_field("Withdrawn", instance.row_count_withdrawn)
        if instance.row_count_delete:
            nb.add_field("Deleted", instance.row_count_delete)
        if instance.row_count_unknown:
            nb.add_field("Unknown", instance.row_count_unknown)
        # and always end with
        if instance.from_file:
            # can only count missing rows if we're from a file and can inject last file data
            nb.add_field("Pre-existing not in import", instance.missing_row_count)

        nb.send()
        if not ongoing_imports:
            classification_imports_complete_signal.send(sender=ClassificationImportRun)
