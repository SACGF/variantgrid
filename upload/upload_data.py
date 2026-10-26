"""
UploadData, the abstract base of the record a processed FileUpload made (UploadedVCF, UploadedGeneList, ...) - one
per FileUpload and file type. Outside upload.models so that an app upload's models import can define its own
(seqauto/models/models_upload.py).
"""
from typing import Optional

from django.db import models

from snpdb.models.models_genome import GenomeBuild


class UploadData(models.Model):
    """ Data created by processing a FileUpload - one per FileUpload/file_type """

    # Retry-import deletes UploadData so the pipeline re-creates it from the file. Set False where the record
    # links to data created before the pipeline started, as the pipeline can't make it again
    created_by_pipeline = True

    class Meta:
        abstract = True

    def get_data(self):
        raise NotImplementedError()

    @property
    def genome_build(self) -> Optional[GenomeBuild]:
        return None

    def get_upload_context(self) -> dict:
        """ Dict for displaying upload widget """
        return {}

    @property
    def requires_user_input(self) -> bool:
        """ The import is waiting for the user to set something the file didn't declare (eg the genome build)
            on the data's page (get_data_url), which finishes the import """
        return False

    def get_data_url(self) -> Optional[str]:
        if data := self.get_data():
            if get_absolute_url := getattr(data, "get_absolute_url", None):
                return get_absolute_url()
        return None
