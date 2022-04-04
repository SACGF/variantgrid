from abc import ABC, abstractmethod
from typing import Iterator
from django.contrib.auth.models import User
from django.db import models
from lazy import lazy
from model_utils.models import TimeStampedModel
from snpdb.models import Lab
import re


class UploadedFileLabStatus(models.TextChoices):
    Pending = 'P', 'Pending'
    Processed = 'MP', 'Processed'
    Error = 'E', 'Error'

    Mapping = "M", "Mapping"
    Importing = "I", "Importing"


# e.g. https://shariant-temp.s3.amazonaws.com/test/hello.txt?AWSAccessKeyId=ASFDWEROMDEOLNZA&Signature=mxBuRkSDFDHgwZyfZYfxQPXE%3D&Expires=1647392831
UPLOADED_FILE_RE = re.compile(r"https:\/\/(?P<bucket>.*?)\.s3\.amazonaws\.com/(?P<file>.*?)(?:\?AWSAccessKeyId.*|$)")


class FileData(ABC):

    @abstractmethod
    def stream(self) -> Iterator[bytes]:
        pass

    @abstractmethod
    def filename(self) -> str:
        pass

    def download_to(self, filename):
        with open(filename, 'wb') as output:
            for chunk in self.stream():
                output.write(chunk)


class FileDataS3(FileData):

    def __init__(self, bucket: str, file: str):
        self.bucket = bucket
        self.file = file

    @property
    def filename(self) -> str:
        segments = self.file.split("/")
        return segments[-1]

    def stream(self) -> Iterator[bytes]:
        from storages.backends.s3boto3 import S3Boto3Storage
        media_storage = S3Boto3Storage(bucket_name=self.bucket)
        with media_storage.open(self.file) as s3file:
            while data := s3file.read(4096):
                yield data


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
    status = models.CharField(max_length=2, choices=UploadedFileLabStatus.choices, default=UploadedFileLabStatus.Pending)

    @lazy
    def file_data(self) -> FileData:
        if match := UPLOADED_FILE_RE.match(self.url):
            bucket = match.group('bucket')
            file = match.group('file')
            return FileDataS3(bucket=bucket, file=file)
        raise ValueError(f"Don't know how to download {self.url}")
