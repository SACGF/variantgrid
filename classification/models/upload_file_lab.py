import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator, List, Dict

import celery
import ijson
from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from lazy import lazy
from model_utils.models import TimeStampedModel

from classification.enums import SubmissionSource
from classification.models.classification_import_run import ClassificationImportRun
from library.utils import batch_iterator
from snpdb.models import Lab
import re
import pathlib


class UploadedFileLabStatus(models.TextChoices):
    Pending = 'P', 'Pending'
    AutoProcessed = 'AP', 'Automatically-Processed'
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


class UploadedFileLab(TimeStampedModel):
    class Meta:
        verbose_name = "Lab Classification File"

    url = models.TextField()
    filename = models.TextField()
    lab = models.ForeignKey(Lab, on_delete=models.CASCADE)
    user = models.ForeignKey(User, on_delete=models.PROTECT)
    comment = models.TextField(default="", blank=True)
    status = models.CharField(max_length=2, choices=UploadedFileLabStatus.choices, default=UploadedFileLabStatus.Pending)

    @lazy
    def file_data(self) -> FileData:
        if match := UPLOADED_FILE_RE.match(self.url):
            bucket = match.group('bucket')
            file = match.group('file')
            return FileDataS3(bucket=bucket, file=file)
        raise ValueError(f"Don't know how to download {self.url}")


class ClassificationImportRunPerformer:

    def __init__(self, upload_file: UploadedFileLab):
        self.upload_file = upload_file

        if omni_importer_dir := settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_APP_DIR:
            self.omni_importer_dir = pathlib.Path(omni_importer_dir)
            if not self.omni_importer_dir.exists():
                raise ValueError(f"{self.omni_importer_dir} does not exist")
            if not (self.omni_importer_dir / "main.py").exists():
                raise ValueError(f"main.py doesn't exist inside OMNI IMPORTER DIR {settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_APP_DIR}")
        else:
            raise ValueError("Require settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_APP_DIR to be set to do automated imports")

        if data_dir := settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_DATA_DIR:
            self.data_dir = pathlib.Path(data_dir)
            self.data_dir.mkdir(parents=True, exist_ok=True)
            if not self.data_dir.exists():
                raise ValueError(f"{self.data_dir} does not exist")
        else:
            raise ValueError(
                "Require settings.VARIANT_CLASSIFICATION_OMNI_IMPORTER_DATA_DIR to be set to do automated imports")

    def update_status(self, status: UploadedFileLabStatus):
        self.upload_file.status = status
        self.upload_file.save()

    @lazy
    def import_id(self) -> str:
        return self.upload_file.file_data.filename + "_" + self.upload_file.lab.group_name

    @celery.shared_task
    def process_async(self):
        self.process()

    @staticmethod
    def cleanup_dir(path: Path, delete_dir: bool = False):
        if existing_files := list(path.iterdir()):
            if len(existing_files) > 4:
                raise ValueError(f"Working subdirectory {path} has {len(existing_files)} files/folders in it, too dangerous to delete")
            if sub_dirs := [ef for ef in existing_files if ef.is_dir()]:
                raise ValueError(
                    f"Working subdirectory {path} has subdirectories {sub_dirs}, too dangerous to delete")

            for sub_file in existing_files:
                sub_file.unlink()
        if delete_dir:
            path.rmdir()

    def process(self):
        from classification.models.classification_inserter import BulkClassificationInserter
        try:
            user = self.upload_file.user
            self.update_status(UploadedFileLabStatus.Mapping)

            file_data = self.upload_file.file_data
            filename = file_data.filename

            sub_folder = f"upload_{self.upload_file.pk}"
            working_sub_folder = self.data_dir / sub_folder
            working_sub_folder.mkdir(exist_ok=True)
            # if working folder already exists, delete files inside it
            ClassificationImportRunPerformer.cleanup_dir(working_sub_folder)

            file_handle = working_sub_folder / filename
            mapped_file = working_sub_folder / f"mapped_classifications.json"

            self.upload_file.file_data.download_to(file_handle)
            # next step is to trigger the omni importer to map the file
            # read the mapped file back
            # and import it
            # --file data/path_west/sharmvl_grch38molecular_variants20220326_060243.json --publish logged_in_users --org path_west --lab unit_1 --env prod
            args: List[str] = [
                settings.PYTHON_COMMAND, "main.py",
                "--file", file_handle,
                "--publish", "logged_in_users",
                "--org", self.upload_file.lab.organization.group_name,
                "--lab", self.upload_file.lab.group_name.split("/")[1],
                "--env", f"file:{mapped_file}"
            ]
            # could make it returned the mapped file to stdout
            # but it's handy having the mapped file present
            process = subprocess.Popen(
                args,
                stdout=subprocess.PIPE,
                cwd=self.omni_importer_dir
            )

            stdout, _ = process.communicate()
            stdout_str = stdout.decode()
            if error_code := process.returncode:
                raise ValueError(f"cwd {self.omni_importer_dir : {' '.join(args)}} returned {error_code}")

            if mapped_file.exists():
                with open(mapped_file, 'r') as file_handle:
                    self.update_status(UploadedFileLabStatus.Importing)

                    def row_generator() -> Dict:
                        nonlocal file_handle
                        for record in ijson.items(file_handle, 'records.item'):
                            yield record

                    for batch in batch_iterator(row_generator(), batch_size=50):
                        # record the import
                        ClassificationImportRun.record_classification_import(
                            identifier=self.import_id,
                            add_row_count=len(batch))

                        bci = BulkClassificationInserter(user=user)
                        for row in batch:
                            bci.insert(data=row, submission_source=SubmissionSource.API)
                        bci.finish()

                    ClassificationImportRun.record_classification_import(
                        identifier=self.import_id,
                        add_row_count=0,
                        is_complete=True
                    )

                    # tidy up if everything is working
                    # if things failed (so we never reached this code), we'll still have this directory to investigate
                    ClassificationImportRunPerformer.cleanup_dir(working_sub_folder, delete_dir=True)

                    self.update_status(UploadedFileLabStatus.Processed)
            else:
                self.update_status(UploadedFileLabStatus.Error)
        except:
            self.update_status(UploadedFileLabStatus.Error)
            raise
