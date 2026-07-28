from dataclasses import dataclass
from datetime import datetime
from typing import Optional, Any

from dateutil import tz
from django.db import models
from django.db.models import CASCADE
from model_utils.models import TimeStampedModel


class ManualMigrationTask(models.Model):
    id = models.TextField(primary_key=True)
    # Gate names (see manual.gates) that must be satisfied before this task may auto-run.
    requires = models.JSONField(null=False, blank=True, default=list)

    @staticmethod
    def describe_manual(text: str) -> str:
        parts = text.split("*", maxsplit=1)
        category = parts[0]
        remainder = parts[1]
        if category == "manage":
            return f"python3 manage.py {remainder}"
        return f"{category} {remainder}"

    def __str__(self):
        return ManualMigrationTask.describe_manual(self.id)


class ManualMigrationRequired(TimeStampedModel):
    task = models.ForeignKey(ManualMigrationTask, on_delete=CASCADE)
    note = models.TextField(null=True, blank=True)


class ManualMigrationAttempt(TimeStampedModel):
    task = models.ForeignKey(ManualMigrationTask, on_delete=CASCADE)
    source_version = models.TextField(null=True, blank=True)
    requires_retry = models.BooleanField(default=False)
    note = models.TextField(null=True, blank=True)


class ManualGateSatisfied(TimeStampedModel):
    """ Records that an operator has confirmed a ManualGate (see manual.gates) as satisfied.
        A gate is considered satisfied while at least one row for its name exists. """
    name = models.TextField()
    source_version = models.TextField(null=True, blank=True)
    note = models.TextField(null=True, blank=True)

    @staticmethod
    def is_satisfied(name: str) -> bool:
        return ManualGateSatisfied.objects.filter(name=name).exists()


@dataclass
class ManualMigrationOutstanding:
    task: ManualMigrationTask
    outstanding_required: list[ManualMigrationRequired]
    last_attempt: Optional[ManualMigrationAttempt]
    last_success: Optional[ManualMigrationAttempt]

    @staticmethod
    def outstanding_task(task: ManualMigrationTask) -> Optional['ManualMigrationOutstanding']:
        last_success = ManualMigrationAttempt.objects.filter(task=task, requires_retry=False).order_by('-created').first()
        last_success_date: datetime
        if last_success:
            last_success_date = last_success.created
        else:
            last_success_date = datetime.min.replace(tzinfo=tz.UTC)

        outstanding_required = ManualMigrationRequired.objects.filter(task=task, created__gte=last_success_date).order_by('created')
        if not outstanding_required:
            return None

        last_attempt = ManualMigrationAttempt.objects.filter(task=task).order_by('-created').first()
        if last_attempt == last_success:
            last_attempt = None  # don't report last attempt twice

        return ManualMigrationOutstanding(
            task=task,
            outstanding_required=list(outstanding_required),
            last_attempt=last_attempt,
            last_success=last_success
        )

    @staticmethod
    def outstanding_tasks() -> list['ManualMigrationOutstanding']:
        outstandings: list[ManualMigrationOutstanding] = []
        for task in ManualMigrationTask.objects.all():
            outstanding = ManualMigrationOutstanding.outstanding_task(task)
            if outstanding:
                outstandings.append(outstanding)
        return outstandings

    def to_json(self) -> dict[str, Any]:
        data: dict[str, Any] = {}
        id_split = self.task.id.split("*", maxsplit=1)
        data["id"] = self.task.id
        data["category"] = id_split[0]
        data["line"] = id_split[1]
        data["requires"] = list(self.task.requires or [])
        data["notes"] = [required.note for required in self.outstanding_required if required.note]
        if self.last_success:
            data["last_success"] = {
                "date": self.last_success.created.timestamp(),
                "note": self.last_success.note
            }
        if self.last_attempt:
            data["last_attempt"] = {
                "date": self.last_attempt.created.timestamp(),
                "note": self.last_attempt.note
            }
        return data
