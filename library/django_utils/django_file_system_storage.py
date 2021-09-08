import os

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.utils.deconstruct import deconstructible


@deconstructible
class PrivateUploadStorage(FileSystemStorage):
    """
        FileSystemStorage hard-codes path in Django migrations - use this to allow relative
        From https://stackoverflow.com/a/32349636/295724 """

    def __init__(self, subdir=''):
        self.subdir = subdir
        super().__init__(location=os.path.join(settings.UPLOAD_DIR, self.subdir), base_url=None)

    def __eq__(self, other):
        return self.subdir == other.subdir
