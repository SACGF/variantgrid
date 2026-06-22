from django.db import models, transaction
from django.utils import timezone
from django_extensions.db.models import TimeStampedModel

__all__ = ["JobsControl"]


class JobsControl(TimeStampedModel):
    """ Singleton operational brake for the analysis + annotation job dispatchers.

        When paused, the single-worker dispatchers (analysis lease_ready_nodes /
        annotation dispatch_annotation_runs) and the stalled-job sweep no-op, so nothing new is
        leased or launched. In-flight node/run statuses are left untouched - resuming lets the
        existing reclaim logic pick work back up where it left off.

        Primary use is a crash safety brake: on a long-lived host a low system uptime means the box
        rebooted (a normal deploy restarts the worker process but leaves uptime high), so the worker
        auto-pauses once per boot (see snpdb.signals.jobs_autopause) - the jobs that may have crashed
        the machine don't immediately re-launch and crash it again. An admin inspects, then resumes
        via the 'jobs_control' management command. """
    SINGLETON_PK = 1

    paused = models.BooleanField(default=False)
    reason = models.TextField(blank=True)
    paused_by = models.TextField(blank=True)
    paused_at = models.DateTimeField(null=True, blank=True)
    # /proc/stat btime (epoch secs) we last auto-paused for, so the reboot auto-pause fires exactly
    # once per boot and doesn't re-trip after an admin resume while uptime is still low.
    last_autopause_boot_time = models.BigIntegerField(null=True, blank=True)

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
    def autopause_for_boot(cls, boot_time: int, reason: str, by: str) -> bool:
        """ Pause once per boot, keyed on the host boot time. Returns True if it paused now, False
            if this boot was already handled (so concurrent workers / restarts don't re-trip it). """
        with transaction.atomic():
            cls.objects.get_or_create(pk=cls.SINGLETON_PK)
            obj = cls.objects.select_for_update().get(pk=cls.SINGLETON_PK)
            if obj.last_autopause_boot_time == boot_time:
                return False
            obj.last_autopause_boot_time = boot_time
            obj.paused = True
            obj.paused_at = timezone.now()
            obj.reason = reason
            obj.paused_by = by
            obj.save()
        return True

    def __str__(self):
        if self.paused:
            return f"JobsControl PAUSED at {self.paused_at} by {self.paused_by}: {self.reason}"
        return "JobsControl: running"
