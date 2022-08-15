import contextlib
import datetime
import re
import tarfile
from abc import ABC, abstractmethod
from io import BytesIO
from os import PathLike
from pathlib import Path
from typing import Iterator, Optional, Union, List
import zipfile
import tarfile
from zipfile import ZipFile


class FileHandle(ABC):

    @abstractmethod
    def _data_handle(self):
        """ Return handle to an open file or storage, anything that you can call .read() on """

    @property
    @abstractmethod
    def filename(self) -> str:
        pass

    @contextlib.contextmanager
    def open(self):
        handle = self._data_handle()
        try:
            yield handle
        finally:
            handle.close()

    def stream(self) -> Iterator[bytes]:
        handle = self._data_handle()
        while buffer := handle.read(4096):
            yield buffer
        handle.close()

    def download_to(self, filename: Union[str, PathLike]):
        with open(filename, 'wb') as output_file:
            with self.open() as input_file:
                output_file.write(input_file.read())

    def download_to_dir(self, download_dir: Path, extract_zip: bool = False):
        with self.open() as input_file:
            if extract_zip and self.filename.endswith(".zip"):
                zippy = ZipFile(input_file)
                zippy.extractall(path=download_dir)
            elif extract_zip and (
                    self.filename.endswith(".tar.gz")
                    or self.filename.endswith(".tar")
                    or self.filename.endswith(".tgz")
            ):
                with tarfile.open(fileobj=input_file, mode="r:*") as tarry:
                    tarry.extractall(path=download_dir)
            else:
                with open(download_dir / self.filename, 'wb') as output_file:
                    output_file.write(input_file.read())

    @abstractmethod
    def save(self, file_obj) -> 'FileHandle':
        pass

    @abstractmethod
    def sub_file(self, sub_file: str) -> 'FileHandle':
        pass

    @abstractmethod
    def list(self) -> List['FileHandle']:
        pass

    @property
    @abstractmethod
    def clean_url(self) -> str:
        pass

    @property
    @abstractmethod
    def modified(self):
        pass

    @property
    @abstractmethod
    def size(self) -> Optional[int]:
        pass

    @property
    @abstractmethod
    def exists(self) -> bool:
        pass


class FileHandleS3(FileHandle):

    def __init__(self, media_storage, file: str):
        self.file = file
        self._media_storage = media_storage

    @staticmethod
    def create(bucket: str, file: str):
        from storages.backends.s3boto3 import S3Boto3Storage
        media_storage = S3Boto3Storage(bucket_name=bucket)
        return FileHandleS3(media_storage=media_storage, file=file)

    def _data_handle(self):
        return self._media_storage.open(self.file)

    @property
    def bucket_name(self):
        return self._media_storage.bucket_name

    @property
    def filename(self) -> str:
        segments = self.file.split("/")
        return segments[-1]

    def sub_path(self, sub_file: str):
        if self.file:
            return self.file + '/' + sub_file
        else:
            return sub_file

    def sub_file(self, sub_file: str) -> 'FileHandleS3':
        return FileHandleS3(media_storage=self._media_storage, file=self.sub_path(sub_file))

    def list(self) -> List['FileHandleS3']:
        _, files = self._media_storage.listdir(self.file)
        return [self.sub_file(file) for file in files]

    def save(self, file_obj) -> 'FileHandleS3':
        file_path = self.sub_path(file_obj.name)
        self._media_storage.save(file_path, file_obj)
        return FileHandleS3(media_storage=self._media_storage, file=file_path)

    @property
    def clean_url(self) -> str:
        # s3 can produce a HTTPS path with temporary access token, but if we want the file path longer term
        # just produce a full s3 file path
        return f"s3://{self.bucket_name}/{self.file}"

    @property
    def modified(self) -> datetime:
        return self._media_storage.get_modified_time(self.file)

    @property
    def size(self) -> Optional[int]:
        return self._media_storage.size(self.file)

    @property
    def exists(self) -> bool:
        return self._media_storage.exists(self.file)

    def __str__(self):
        return f"{self.clean_url} {self.modified} {self.size/1024}KB"


UPLOADED_S3_TEMP_URL = re.compile(r"https:\/\/(?P<bucket>.*?)\.s3\.amazonaws\.com/(?P<file>.*?)(?:\?AWSAccessKeyId.*|$)")
UPLOADED_S3_CLEAN_URL = re.compile(r"s3:\/\/(?P<bucket>.*?)\/(?P<file>.*)")


def resolve_uploaded_url_to_handle(url: str) -> FileHandle:
    # e.g. https://shariant-temp.s3.amazonaws.com/test/hello.txt?AWSAccessKeyId=ASFDWEROMDEOLNZA&Signature=mxBuRkSDFDHgwZyfZYfxQPXE%3D&Expires=1647392831
    if match := UPLOADED_S3_TEMP_URL.match(url) or UPLOADED_S3_CLEAN_URL.match(url):
        bucket = match.group('bucket')
        file = match.group('file')
        return FileHandleS3.create(bucket=bucket, file=file)
    raise ValueError(f"Could not parse URL - or not of supported type {url}")
