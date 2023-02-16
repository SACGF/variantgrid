from functools import cached_property
from typing import Union, Optional, Dict, List, Tuple, Any, Iterable

from django.contrib.auth.models import User
from django.db import models
from django.urls import reverse
from model_utils.models import TimeStampedModel

from classification.models.uploaded_file_types import FileHandle, resolve_uploaded_url_to_handle
from library.utils import ExportRow, export_column
from snpdb.models import Lab


class UploadedClassificationsUnmappedStatus(models.TextChoices):
    Manual = "MA", "Manual Review Pending"
    Pending = 'P', 'Pending'

    Mapping = "M", "Mapping"
    Importing = "I", "Importing"

    Validated = 'V', 'Validated'
    Error = 'E', 'Error'
    Processed = 'MP', 'Processed'


class UploadedClassificationsUnmappedValidationRow(ExportRow):

    def __init__(self, row: Dict[str, Any]):
        self._row = row

    def filename(self):
        return self._row.get('file')

    @export_column("line number")
    def line_number(self) -> Optional[int]:
        try:
            return int(self._row.get('row'))
        except:
            return None

    @export_column("severity")
    def severity(self):
        return self._row.get('severity')

    @property
    def filename_line_number(self) -> str:
        return f"{self.filename() or '?'}:{self.line_number() or '?'}"

    @export_column("category")
    def category(self):
        return self._row.get('category')

    @export_column("message")
    def message(self):
        return self._row.get('message')


class UploadedClassificationsUnmapped(TimeStampedModel):
    class Meta:
        verbose_name = "Classification upload file"

    url = models.TextField()
    filename = models.TextField()
    lab = models.ForeignKey(Lab, on_delete=models.CASCADE)
    user = models.ForeignKey(User, on_delete=models.PROTECT)
    comment = models.TextField(default="", blank=True)
    validation_summary = models.JSONField(null=True, blank=True)
    validation_list = models.JSONField(null=True, blank=True)
    status = models.CharField(max_length=2, choices=UploadedClassificationsUnmappedStatus.choices, default=UploadedClassificationsUnmappedStatus.Pending)

    # taken from the server details
    effective_modified = models.DateTimeField(null=True, blank=True)
    file_size = models.IntegerField(null=True, blank=True)

    def __str__(self):
        return f"(file_id={self.pk}) {self.lab} {self.filename} {self.created}"

    def get_absolute_url(self):
        return reverse('classification_upload_unmapped_status', kwargs={'uploaded_classification_unmapped_id': self.pk})

    def validation_list_objs(self) -> List[UploadedClassificationsUnmappedValidationRow]:
        if messages := self.validation_list.get('messages'):
            return [UploadedClassificationsUnmappedValidationRow(entry) for entry in messages]
        return []

    @property
    def validation_includes_json(self) -> bool:
        if self.filename.endswith(".json"):
            return True
        elif self.filename.endswith(".zip"):
            try:
                if messages := self.validation_list.get('messages'):
                    for message in messages:
                        if message.get('file').endswith('.json'):
                            return True
            except Exception:
                pass
        return False

    @cached_property
    def file_data(self) -> FileHandle:
        return resolve_uploaded_url_to_handle(self.url)

    @property
    def message_counts(self) -> Optional[List[Tuple[str, int]]]:
        if summary := self.validation_summary:
            entries = list(summary.get("message_counts").items())
            return sorted(entries, key=lambda x: x[0])

    @property
    def validation_summary_properties(self) -> Optional[List[Tuple[str, Union[int, str]]]]:
        if summary := self.validation_summary:
            entries = [(key, value) for key, value in summary.items() if isinstance(value, (str, int))]
            return sorted(entries, key=lambda x: x[0])
