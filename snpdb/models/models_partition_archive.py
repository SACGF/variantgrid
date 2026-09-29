from django.contrib.auth.models import User
from django.contrib.postgres.fields import ArrayField
from django.db import models
from django.db.models.deletion import SET_NULL


class PartitionArchive(models.Model):
    """ Tracks a pg_dump+drop of a RelatedModelsPartitionModel's child partition tables.

        @see claude/issue_1537_archive_plan.md
    """

    class Status(models.TextChoices):
        PENDING = "PENDING", "Pending"
        IN_PROGRESS = "IN_PROGRESS", "In progress"
        COMPLETE = "COMPLETE", "Complete"
        FAILED = "FAILED", "Failed"

    archive_name = models.TextField(unique=True)
    dump_path = models.TextField()
    sha256 = models.TextField(null=True, blank=True)

    source_app_label = models.TextField()
    source_model = models.TextField()
    source_pk = models.TextField()
    source_table_names = ArrayField(models.TextField())

    status = models.CharField(max_length=20, choices=Status.choices, default=Status.PENDING)
    error_message = models.TextField(null=True, blank=True)

    dumped_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL,
                                  related_name="dumped_partition_archives")
    dumped_date = models.DateTimeField(auto_now_add=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)

    cleared_date = models.DateTimeField(null=True, blank=True)
    cleared_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL,
                                   related_name="cleared_partition_archives")

    restored_date = models.DateTimeField(null=True, blank=True)
    restored_by = models.ForeignKey(User, null=True, blank=True, on_delete=SET_NULL,
                                    related_name="restored_partition_archives")

    class Meta:
        ordering = ("-dumped_date",)

    def __str__(self):
        return f"{self.archive_name} ({self.get_status_display()})"

    @property
    def source_label(self) -> str:
        return f"{self.source_app_label}.{self.source_model}#{self.source_pk}"

    def resolve_source(self):
        """ Return the live source row if it still exists, else None. """
        from django.apps import apps
        model = apps.get_model(self.source_app_label, self.source_model)
        return model.objects.filter(pk=self.source_pk).first()
