from typing import NamedTuple, Optional

from django.db import models, transaction
from django.utils import timezone
from django_extensions.db.models import TimeStampedModel

__all__ = ["JobsControl", "BootRegistration"]


class BootRegistration(NamedTuple):
    """ new_boot is False when this host boot was already recorded (concurrent workers / later
        worker restarts within the one boot). previous_boot_time is the boot before this one, or
        None if this is the first boot we've recorded. """
    new_boot: bool
    previous_boot_time: Optional[int]


class JobsControl(TimeStampedModel):
    """ Singleton operational brake for the analysis + annotation job dispatchers.

        When paused, the single-worker dispatchers (analysis lease_ready_nodes /
        annotation dispatch_annotation_runs) and the stalled-job sweep no-op, so nothing new is
        leased or launched. In-flight node/run statuses are left untouched - resuming lets the
        existing reclaim logic pick work back up where it left off.

        Primary use is a crash safety brake: the worker records each host boot it starts under, and
        auto-pauses when two boots land close together (see snpdb.signals.jobs_autopause) - a box
        that reboots, comes up and then goes down again is likely being crashed by the jobs it
        relaunches, so we hold them until an admin inspects and resumes via the 'jobs_control'
        management command. A single reboot the box survives keeps running. """
    SINGLETON_PK = 1

    paused = models.BooleanField(default=False)
    reason = models.TextField(blank=True)
    paused_by = models.TextField(blank=True)
    paused_at = models.DateTimeField(null=True, blank=True)
    # The two most recent /proc/stat btime values (epoch secs) a worker has started under. The gap
    # between them is what tells a one-off reboot from a crash loop, and keying on the latest means
    # each boot is judged exactly once - later worker restarts don't re-trip after an admin resume.
    last_boot_time = models.BigIntegerField(null=True, blank=True)
    previous_boot_time = models.BigIntegerField(null=True, blank=True)

    def save(self, *args, **kwargs):
        self.pk = self.SINGLETON_PK  # enforce a single row
        super().save(*args, **kwargs)

    @classmethod
    def get_solo(cls) -> "JobsControl":
        obj, _ = cls.objects.get_or_create(pk=cls.SINGLETON_PK)
        return obj

    @classmethod
    def is_paused(cls) -> bool:
        """ Cheap indexed PK check - called on every dispatch. Missing row => not paused. """
        return cls.objects.filter(pk=cls.SINGLETON_PK, paused=True).exists()

    @classmethod
    def pause(cls, reason: str = "", by: str = "") -> "JobsControl":
        with transaction.atomic():
            obj = cls.get_solo()
            if not obj.paused:
                obj.paused_at = timezone.now()
            obj.paused = True
            obj.reason = reason
            obj.paused_by = by
            obj.save()
        return obj

    @classmethod
    def resume(cls, by: str = "") -> "JobsControl":
        with transaction.atomic():
            obj = cls.get_solo()
            obj.paused = False
            obj.reason = ""
            obj.paused_by = by
            obj.paused_at = None
            obj.save()
        return obj

    @classmethod
    def register_boot(cls, boot_time: int) -> BootRegistration:
        """ Record the host boot this worker started under, shuffling the previous one down. """
        with transaction.atomic():
            cls.objects.get_or_create(pk=cls.SINGLETON_PK)
            obj = cls.objects.select_for_update().get(pk=cls.SINGLETON_PK)
            if obj.last_boot_time == boot_time:
                return BootRegistration(new_boot=False, previous_boot_time=obj.previous_boot_time)
            previous_boot_time = obj.last_boot_time
            obj.previous_boot_time = previous_boot_time
            obj.last_boot_time = boot_time
            obj.save()
        return BootRegistration(new_boot=True, previous_boot_time=previous_boot_time)

    def __str__(self):
        if self.paused:
            return f"JobsControl PAUSED at {self.paused_at} by {self.paused_by}: {self.reason}"
        return "JobsControl: running"
